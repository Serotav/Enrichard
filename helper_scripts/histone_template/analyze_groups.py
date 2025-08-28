#!/usr/bin/env python3

import argparse
import pathlib
import polars as pl
import numpy as np
import sys
from scipy.stats import ttest_ind

def _aggregate_and_melt_heatmaps(directory: pathlib.Path) -> pl.DataFrame:
    """
    Reads all wide-format TSV heatmap files in a directory, melts them
    into a long format, and concatenates them into a single DataFrame using
    optimized, batch-processing methods.
    """
    all_files = list(directory.glob('*.tsv'))
    if not all_files:
        return pl.DataFrame()

    # Lazily scan all files at once. This is much faster than a Python loop.
    # We will melt them individually within a list comprehension and then concat.
    lazy_frames = []
    for file_path in all_files:
        try:
            # Lazily read each file to get its schema for the melt operation.
            lf = pl.scan_csv(file_path, separator='\t', infer_schema_length=1)
            
            # The first column name is the identifier for the melt operation.
            id_col_name = lf.columns[0]
            
            # Create a lazy plan for melting this specific file.
            melted_lf = lf.melt(
                id_vars=[id_col_name],
                variable_name="State",
                value_name="Odds-Ratio"
            ).rename({id_col_name: "Cell_Type"})
            
            lazy_frames.append(melted_lf)
        except Exception as e:
            print(f"Warning: Skipping file {file_path} due to an error: {e}")
            continue
            
    if not lazy_frames:
        return pl.DataFrame()

    # Concatenate all the lazy plans and execute them together.
    return pl.concat(lazy_frames).collect()



def run_ttest_analysis(real_df: pl.DataFrame, control_df: pl.DataFrame, n_real_total: int, n_control_total: int, p_cutoff: float) -> pl.DataFrame:
    """
    Performs group-level enrichment on heatmap data using a t-test with
    a robust, optimized Polars approach.
    """
    print("Running group comparison using vectorized T-test...")

    # Step 1: Combine real and control data into a single DataFrame.
    combined_df = pl.concat([
        real_df.with_columns(pl.lit("real").alias("group")),
        control_df.with_columns(pl.lit("control").alias("group"))
    ])

    # Step 2: Group by the composite key ("Cell_Type", "State") and aggregate ORs into lists.
    analysis_df = combined_df.group_by("Cell_Type", "State").agg(
        pl.col("Odds-Ratio").filter(pl.col("group") == "real").alias("real_ors"),
        pl.col("Odds-Ratio").filter(pl.col("group") == "control").alias("control_ors")
    )

    # Step 3: Define a robust function to perform the t-test on the lists.
    def perform_ttest_on_lists(row):
        real_ors = row["real_ors"] if row["real_ors"] is not None else []
        control_ors = row["control_ors"] if row["control_ors"] is not None else []
        
        real_ors.extend([1.0] * (n_real_total - len(real_ors)))
        control_ors.extend([1.0] * (n_control_total - len(control_ors)))
        
        if len(real_ors) < 2 or len(control_ors) < 2:
            return {"T-statistic": np.nan, "P-value": np.nan, "Mean_logOR_Real": np.nan, "Mean_logOR_Control": np.nan}

        log_real_ors = np.log(real_ors)
        log_control_ors = np.log(control_ors)
        stat, p_value = ttest_ind(log_real_ors, log_control_ors, equal_var=False, nan_policy='omit')
        
        return {
            "T-statistic": stat, 
            "P-value": p_value, 
            "Mean_logOR_Real": np.mean(log_real_ors), 
            "Mean_logOR_Control": np.mean(log_control_ors)
        }

    # Step 4: Apply the function across all groups at once.
    results_df = analysis_df.with_columns(
        pl.struct(["real_ors", "control_ors"]).map_elements(
            perform_ttest_on_lists,
            return_dtype=pl.Struct([
                pl.Field("T-statistic", pl.Float64), 
                pl.Field("P-value", pl.Float64),
                pl.Field("Mean_logOR_Real", pl.Float64),
                pl.Field("Mean_logOR_Control", pl.Float64)
            ])
        ).alias("ttest_result")
    ).unnest("ttest_result").drop(["real_ors", "control_ors"]).drop_nulls()

    # Step 5: Filter by the p-value cutoff.
    return results_df.filter(pl.col("P-value") < p_cutoff)


def main():
    """
    Main function to run the group comparison analysis for chromatin heatmaps.
    """
    parser = argparse.ArgumentParser(description="Perform group-level enrichment analysis on chromatin heatmap data using t-tests.")
    parser.add_argument("--analysis-dir", type=pathlib.Path, required=True, help="Path to the module's analysis directory containing 'real_results' and 'control_results' subfolders.")
    parser.add_argument("--p-value-cutoff", type=float, default=0.05, help="P-value cutoff to filter the final group-level results.")
    args = parser.parse_args()

    real_dir = args.analysis_dir / "real_results"
    control_dir = args.analysis_dir / "control_results"
    output_dir = args.analysis_dir
    
    output_filename = "group_comparison_heatmap_ttest.tsv"
    output_path = output_dir / output_filename
    
    n_real_total = len(list(real_dir.glob('*.tsv')))
    n_control_total = len(list(control_dir.glob('*.tsv')))

    if n_real_total == 0 or n_control_total == 0:
        print("Error: 'real_results' and 'control_results' directories must both contain at least one file.", file=sys.stderr)
        return

    # Aggregate and transform data from both groups using the optimized function
    print("Aggregating results from real samples...")
    real_long_df = _aggregate_and_melt_heatmaps(real_dir)
    print("Aggregating results from control samples...")
    control_long_df = _aggregate_and_melt_heatmaps(control_dir)

    if real_long_df.is_empty() and control_long_df.is_empty():
        print("Warning: No valid data found in either real or control directories.")
        return

    # Run the fully optimized t-test analysis, replacing the old loop
    final_df = run_ttest_analysis(real_long_df, control_long_df, n_real_total, n_control_total, args.p_value_cutoff)

    if final_df.is_empty():
        print(f"Warning: No (Cell_Type, State) pairs passed the significance threshold of p < {args.p_value_cutoff}.")
        output_path.touch()
        return

    # Save the final long-format results
    final_df.sort("Cell_Type").write_csv(output_path, separator='\t')
    print(f"Successfully wrote {final_df.height} significant results to {output_path}")

if __name__ == "__main__":
    main()