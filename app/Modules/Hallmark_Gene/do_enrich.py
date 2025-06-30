import polars as pl
from scipy.stats import fisher_exact
import sys
import argparse
import os
import pathlib
from pathlib import Path
from statsmodels.stats.multitest import multipletests

APP_DIR = pathlib.Path(__file__).parent
BACKGROUND_DIR  = APP_DIR/"background"
MERGE_PATH = BACKGROUND_DIR / "merge.parquet" 
PROBE_ID_COL = "Probe_ID"

USER_CUSTOM_BACKGROUND_NAME = os.getenv("USER_CUSTOM_BACKGROUND_NAME")
USER_SAMPLE_NAME = os.getenv("USER_SAMPLE_NAME")
OUTPUT_FILE_NAME = str(APP_DIR.name) + ".tsv"

BITS_PER_CHUNK = 64


def perform_enrichment_from_cache(cache_df: pl.DataFrame,background_df: pl.DataFrame, sample_probes_df: pl.DataFrame) -> pl.DataFrame:
    """
    Performs Fisher's exact test using a cached traits counts.
    """
    probe_id_dtype = background_df.get_column(PROBE_ID_COL).dtype
    sample_ids = sample_probes_df.select(
        pl.col(PROBE_ID_COL).unique().cast(probe_id_dtype),
    )

    sample_df = background_df.join(
        sample_ids, on=PROBE_ID_COL, how="inner"
    )

    contingency_table_df = sample_df.select(
        [pl.col(col).sum().alias(col)for col in sample_df.columns if col != PROBE_ID_COL]
    ).unpivot(index=[], variable_name="Trait", value_name="a")

    total_sample_size = sample_df.height
    total_background_only_size = background_df.height - total_sample_size

    results_df = cache_df.lazy().join(contingency_table_df.lazy(), on='Trait', how="inner").with_columns(
        a = pl.col('a'),
        c = pl.col('totals') - pl.col('a')
    ).drop('totals').filter(
        (pl.col("a") + pl.col("c")) > 0
    ).with_columns(
        b=pl.lit(total_sample_size) - pl.col("a"),
        d=pl.lit(total_background_only_size) - pl.col("c")
    ).collect()

    if results_df.is_empty():
        return pl.DataFrame()


    if results_df.is_empty():
        return pl.DataFrame()

    # Run Fisher's Test and Unnest the results
    def run_fisher(s):
        odds_ratio, p_value = fisher_exact([[s["a"], s["b"]], [s["c"], s["d"]]], alternative='greater')
        return p_value

    freq_in_sample = pl.col("a") / total_sample_size
    freq_in_background = pl.col("c") / total_background_only_size

    final_df = results_df.with_columns(
        pl.struct(["a", "b", "c", "d"]).map_elements(
            run_fisher,
            return_dtype=pl.Float64
        ).alias("P-Value"),

        pl.when(freq_in_background > 0)
          .then(freq_in_sample / freq_in_background)
          .otherwise(None) # Set to null if background frequency is 0 to avoid 'inf'
          .alias("Fold-Change")
    
    )

    return final_df

def perform_enrichment(background_df: pl.DataFrame, sample_probes_df: pl.DataFrame) -> pl.DataFrame:
    """
    Performs Fisher's exact test using a fully vectorized Polars backend.
    """
   
    # Crete the sample df, and the background df with the `is_in_sample` column
    probe_id_dtype = background_df.get_column(PROBE_ID_COL).dtype
    sample_ids = sample_probes_df.select(
        pl.col(PROBE_ID_COL).unique().cast(probe_id_dtype),
        pl.lit(True).alias("is_in_sample_marker")
    )

    background_df = background_df.join(
        sample_ids, on=PROBE_ID_COL, how="left"
    ).with_columns(
        is_in_sample=pl.col("is_in_sample_marker").fill_null(False)
    ).drop("is_in_sample_marker")


    total_sample_size = background_df.filter(pl.col("is_in_sample")).height
    total_background_only_size = background_df.height - total_sample_size
    
    if total_sample_size == 0:
        print("Warning: No probes from the sample were found in the background set.", file=sys.stderr)
        return pl.DataFrame()

    contingency_table_df = background_df.group_by("is_in_sample"
    ).agg([pl.col(col).sum().alias(col) for col in background_df.columns if col != PROBE_ID_COL and col != "is_in_sample"]
    ).unpivot(index="is_in_sample", variable_name="Trait", value_name="count"
    ).pivot(index="Trait", on="is_in_sample", values="count"
    ).rename({"true": "a", "false": "c"}
    ).fill_null(0)
    
    
    results_df = contingency_table_df.filter(
        (pl.col("a") + pl.col("c")) > 0
    ).with_columns(
        b=pl.lit(total_sample_size) - pl.col("a"),
        d=pl.lit(total_background_only_size) - pl.col("c")
    )
    
    results_df = contingency_table_df.filter(
        (pl.col("a") + pl.col("c")) > 0
    ).with_columns(
        b=pl.lit(total_sample_size) - pl.col("a"),
        d=pl.lit(total_background_only_size) - pl.col("c")
    )
    
    if results_df.is_empty():
        return pl.DataFrame()

    # Run Fisher's Test and Unnest the results
    def run_fisher(s):
        odds_ratio, p_value = fisher_exact([[s["a"], s["b"]], [s["c"], s["d"]]], alternative='greater')
        return p_value

    freq_in_sample = pl.col("a") / total_sample_size
    freq_in_background = pl.col("c") / total_background_only_size

    final_df = results_df.with_columns(
        pl.struct(["a", "b", "c", "d"]).map_elements(
            run_fisher,
            return_dtype=pl.Float64
        ).alias("P-Value"),

        pl.when(freq_in_background > 0)
          .then(freq_in_sample / freq_in_background)
          .otherwise(None) # Set to null if background frequency is 0 to avoid 'inf'
          .alias("Fold-Change")
    
    )

    print("Final DataFrame after Fisher's test:")
    return final_df

def apply_correction(results_df: pl.DataFrame, method: str, p_value_threshold: float) -> pl.DataFrame:
    """
    Applies a multiple testing correction to a DataFrame of enrichment results.

    Args:
        results_df (pl.DataFrame): The raw, uncorrected results from perform_enrichment.
        method (str): The correction method to use ('fdr_bh', 'bonferroni', or 'none').
        p_value_threshold (float): The significance threshold.

    Returns:
        pl.DataFrame: A filtered DataFrame containing only the significant results after correction.
    """
    if results_df.is_empty():
        return results_df

    if method == 'none':
        # Filter by the original P-Value and create a new column for P-adj
        return results_df.filter(
            pl.col("P-Value") < p_value_threshold
        ).with_columns(
            pl.col("P-Value").alias("P-adj")
        )

    # the multipletests function only works with a list of p-values, so we extract them
    p_values = results_df["P-Value"]
    
    # Use statsmodels to perform the correction
    try:
        reject, pvals_corrected, _, _ = multipletests(
            pvals=p_values, 
            alpha=p_value_threshold, 
            method=method, 
            is_sorted=False
        )   
    except ValueError as e:
        print(f"Error applying multiple testing correction: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Add the adjusted p-values as a new column
    corrected_df = results_df.with_columns(
        pl.Series(name="P-adj", values=pvals_corrected)
    )
    
    # Filter the DataFrame based on the new adjusted p-values
    significant_results = corrected_df.filter(
        pl.col("P-adj") < p_value_threshold
    )
    
    return significant_results.sort("P-adj")

def merge_background():
    
    source_files = [
        file_path for file_path in BACKGROUND_DIR.glob("*.parquet")
        if file_path.name != MERGE_PATH.name # Exclude the merge file itself, WE SHOULD NOT BE HERE IF IT EXISTS
    ]

    if not source_files:
        print(f"No files to merge in {BACKGROUND_DIR} directory.", file=sys.stderr)
        return pl.DataFrame() # Return empty DataFrame if no files

    lazy_frames = [
        pl.scan_parquet(file_path)
        for file_path in source_files
    ]

    #  Concatenate all lazy frames, drop duplicates, and then collect the result.
    merged_df = pl.concat(lazy_frames).unique(subset=[PROBE_ID_COL], keep='first').collect()

    merged_df.write_parquet(MERGE_PATH)
    return merged_df

def custom_background(user_background_file: Path) -> pl.DataFrame:
    if not MERGE_PATH.exists():
        merge_df = merge_background()
    else: 
        merge_df = pl.read_parquet(MERGE_PATH)
    
    user_probes_df = pl.read_csv(user_background_file, has_header=False).rename({"column_1": PROBE_ID_COL})
    
    custom_annotated_df = merge_df.join(user_probes_df, on=PROBE_ID_COL, how="inner")
    
    return custom_annotated_df

def get_background_df(background_name: str, user_dir: Path) -> pl.DataFrame:
    if background_name == USER_CUSTOM_BACKGROUND_NAME:
        print("Processing custom background...", file=sys.stderr)
        # The user's custom background file was saved in their session directory
        user_background_file = user_dir / USER_CUSTOM_BACKGROUND_NAME
        return custom_background(user_background_file)
    
    # Find a file in BACKGROUND_DIR that starts with args.background_name
    background_name += "."
    background_dir = pathlib.Path(BACKGROUND_DIR)
    matching_files = list(background_dir.glob(f"{background_name}*"))
    if not matching_files or len(matching_files) > 1:
        print(f"Error: No file found in {background_dir} starting with '{background_name}' or multiple found files found:[{len(matching_files)}]", file=sys.stderr)
        sys.exit(1)
    
    # Read with the right types
    file_path = matching_files[0]

    return pl.read_parquet(file_path)


def main()-> None:

    parser = argparse.ArgumentParser(description="Perform Fisher's exact test for enrichment analysis on CpG sites.")
    parser.add_argument("--user_dir", help="Path to the user directory.")
    parser.add_argument("--background_name", help="Name of the background annotation file.")
    parser.add_argument("--p_value", type=float, default=0.05, help="P-value threshold for significance.")
    parser.add_argument("--correction", default='none',  help="Multiple testing correction method to apply.")
    parser.add_argument("--output_folder", help="Path to the output folder where results will be saved.")
    parser.add_argument("--cache_folder", help="Path to the cache folder.")

    args = parser.parse_args()
    user_dir = Path(args.user_dir)
    print(f"User directory: {user_dir}", file=sys.stderr)
    # Load Data 
    try:
        background_df = get_background_df(args.background_name, user_dir)
        sample_df = pl.read_csv(Path(args.user_dir) / USER_SAMPLE_NAME, has_header=False).rename({ "column_1": PROBE_ID_COL })
    except Exception as e:
        print(f"Error reading background or sample {e}", file=sys.stderr)
        sys.exit(1)

    # Run Analysis
    if args.background_name== USER_CUSTOM_BACKGROUND_NAME:
        raw_results_df = perform_enrichment(background_df, sample_df)
    else:
        cache_path = list(Path(args.cache_folder).glob(f"{args.background_name}.*"))
        if not cache_path or len(cache_path) > 1:
            print(f"Error: No cache file found for background '{args.background_name}' or multiple found [{len(cache_path)}].", file=sys.stderr)
            sys.exit(1)
        cache_df = pl.read_parquet(cache_path[0])
        raw_results_df = perform_enrichment_from_cache(cache_df, background_df, sample_df)


    final_results_df = apply_correction(raw_results_df, args.correction, args.p_value)
    # Save Results 
    # Check if output directory exists, create if not
    output_dir = pathlib.Path(args.output_folder) / APP_DIR.name
    output_dir.mkdir(parents=True, exist_ok=True)
    destination = output_dir/ OUTPUT_FILE_NAME
    
    final_results_df.write_csv(destination, separator="\t")


if __name__ == "__main__":
    main()