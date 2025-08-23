import polars as pl
from scipy.stats import fisher_exact
import sys
import argparse
import os
import pathlib
from pathlib import Path
from statsmodels.stats.multitest import multipletests

APP_DIR = pathlib.Path(__file__).parent
PROBE_ID_COL = "Probe_ID"
USER_SAMPLE_NAME = os.getenv("USER_SAMPLE_NAME")

def get_corrected_p_values(results_df: pl.DataFrame, method: str, p_value_threshold: float) -> pl.DataFrame:
    """Applies multiple testing correction and returns the full dataframe with an added P-adj column."""
    if results_df.is_empty() or method == 'none':
        return results_df.with_columns(pl.col("P-Value").alias("P-adj"))

    p_values = results_df["P-Value"].to_numpy()
    try:
        _reject, pvals_corrected, _, _ = multipletests(pvals=p_values, alpha=p_value_threshold, method=method)
        return results_df.with_columns(
            pl.Series(name="P-adj", values=pvals_corrected)
        )
    except Exception as e:
        print(f"Error during correction: {e}", file=sys.stderr)
        # Return with a null P-adj column on error
        return results_df.with_columns(pl.lit(None, dtype=pl.Float64).alias("P-adj"))

def run_fisher_test(contingency_df: pl.DataFrame) -> pl.DataFrame:
    """Runs Fisher's exact test on a DataFrame containing contingency table values."""
    if contingency_df.is_empty():
        return pl.DataFrame()

    def run_fisher_struct(row):
        odds, pval = fisher_exact([[row["a"], row["b"]], [row["c"], row["d"]]], alternative='greater')
        return {"Odds-Ratio": odds, "P-Value": pval}

    results_struct_col = pl.struct(["a", "b", "c", "d"]).map_elements(
        run_fisher_struct,
        return_dtype=pl.Struct([pl.Field("Odds-Ratio", pl.Float64), pl.Field("P-Value", pl.Float64)])
    ).alias("results")

    return contingency_df.with_columns(results_struct_col).unnest("results")

def main():
    parser = argparse.ArgumentParser(description="Perform grouped chromatin state enrichment analysis.")
    parser.add_argument("--user_dir", required=True, help="Path to the user directory.")
    parser.add_argument("--background_name", required=True, help="Name of the background set (e.g., 'HM450').")
    parser.add_argument("--p_value", type=float, default=0.05, help="P-value threshold.")
    parser.add_argument("--correction", default='none', help="Multiple testing correction method.")
    parser.add_argument("--output_folder", required=True, help="Path to the output folder.")
    args = parser.parse_args()

    # Load Data
    try:
        background_path = next((APP_DIR / "background").glob(f"{args.background_name}.hg*annotated.parquet"))
        background_df = pl.read_parquet(background_path)
        sample_probes_df = pl.read_csv(Path(args.user_dir) / USER_SAMPLE_NAME, has_header=False).rename({"column_1": PROBE_ID_COL})
        print(f"chromatine module, using {background_path} for {args.user_dir}")
    except (StopIteration, Exception) as e:
        print(f"Error loading data files for {args.background_name}: {e}", file=sys.stderr)
        sys.exit(1)

    # Prepare for Analysis
    all_traits = [col for col in background_df.columns if col != PROBE_ID_COL]
    
    sample_annotated_df = background_df.join(sample_probes_df, on=PROBE_ID_COL, how="inner")
    total_sample_size = len(sample_annotated_df)
    
    if total_sample_size == 0:
        print("Warning: No probes from your sample were found in the background.", file=sys.stderr)
        sys.exit(0)
        
    total_background_size = len(background_df)
    total_background_only_size = total_background_size - total_sample_size

    # Calculate Contingency Table
    print("Calculating contingency table...", file=sys.stderr)
    
    sample_counts = sample_annotated_df.select(all_traits).sum().unpivot(
        index=[], variable_name="Trait", value_name="a"
    )
    
    total_counts = background_df.select(all_traits).sum().unpivot(
        index=[], variable_name="Trait", value_name="a_plus_c"
    )

    contingency_df = sample_counts.join(total_counts, on="Trait", how="left").with_columns(
        b = total_sample_size - pl.col('a'),
        c = pl.col('a_plus_c') - pl.col('a'),
        d = total_background_only_size - (pl.col('a_plus_c') - pl.col('a'))
    ).drop("a_plus_c")

    #  Run Analysis 
    print("Running enrichment tests...", file=sys.stderr)
    raw_results_df = run_fisher_test(contingency_df)
    
    print("Applying multiple testing correction...", file=sys.stderr)
    corrected_results_df = get_corrected_p_values(raw_results_df, args.correction, args.p_value)

    #  Reshape for Heatmap Output 
    output_dir = Path(args.output_folder) / APP_DIR.name
    output_dir.mkdir(parents=True, exist_ok=True)
    destination = output_dir / f"{APP_DIR.name}_HEATMAP.tsv"

    if corrected_results_df.is_empty():
        print("No results to process.", file=sys.stderr)
        destination.touch()
        sys.exit(0)

    print("Reshaping results into a heatmap...", file=sys.stderr)
    
    # Set Odds-Ratio to 1.0 for non-significant results, then pivot
    heatmap_df = corrected_results_df.with_columns(
        Display_OR = pl.when(pl.col("P-adj") < args.p_value)
                       .then(pl.col("Odds-Ratio"))
                       .otherwise(1.0)
    ).with_columns(
        pl.col("Trait").str.split_exact("_", 1).struct.rename_fields(["Cell_Type", "State"])
    ).unnest("Trait")

    pivoted_df = heatmap_df.pivot(
        index="Cell_Type",
        columns="State",
        values="Display_OR"
    ).fill_null(1.0)

    # Add human-readable names and save
    print("Mapping EID to cell type names...", file=sys.stderr)
    try:
        script_dir = pathlib.Path(__file__).parent
        metadata_path = script_dir / "EID_metadata.tab"
        metadata_df = pl.read_csv(metadata_path, separator='\t').select(["EID", "STD_NAME"])

        # Get the state columns before the join
        state_columns = pivoted_df.columns[1:]

        final_df = pivoted_df.join(
            metadata_df, left_on="Cell_Type", right_on="EID", how="left"
        ).select(
            pl.col("STD_NAME").alias("Cell_Type"), # Use the human-readable name
            *state_columns # Keep all the state columns
        ).fill_null("Unknown")

        final_df.write_csv(destination, separator="\t")
        print(f"Chromatine state analysis complete. Heatmap saved to {destination}", file=sys.stderr)

    except Exception as e:
        print(f"Failed to map EIDs to names. Saving with EIDs. Error: {e}", file=sys.stderr)
        pivoted_df.write_csv(destination, separator="\t")


if __name__ == "__main__":
    main()
