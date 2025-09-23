#!/usr/bin/env python3

import argparse
import pathlib
import polars as pl
import numpy as np
from scipy.stats import ttest_ind

def _aggregate_results(directory: pathlib.Path) -> pl.DataFrame:
    """
    Finds and concatenates all result .tsv files in a given directory using
    Polars' lazy evaluation for optimized batch processing. It calculates
    the Odds-Ratio and selects the required columns.
    """
    # Define the glob pattern to find all .tsv files in the directory.
    file_pattern = str(directory / '*.tsv')

    required_schema = {
        "Trait": pl.String,
        "a": pl.Int64,
        "b": pl.Int64,
        "c": pl.Int64,
        "d": pl.Int64,
        "P-adj": pl.Float64,
        "P-Value": pl.Float64,
        "Fold-Change": pl.Float64
    }
    
    try:
        lazy_df = pl.scan_csv(
            file_pattern, 
            separator='\t', 
            null_values="NA", 
            schema=required_schema,
            infer_schema_length=0 
        )

        processed_lazy_df = (
            lazy_df
            .rename({"P-adj": "P-value"})
            .with_columns(
                Odds_Ratio=((pl.col("a") + 0.5) * (pl.col("d") + 0.5)) / ((pl.col("b") + 0.5) * (pl.col("c") + 0.5))
            )
            .rename({"Odds_Ratio":"Odds-Ratio"})
            .select(["Trait", "Odds-Ratio", "P-value"])
        )

        return processed_lazy_df.collect()

    except Exception as e:
        print(f"An error occurred during batch processing of files in {directory}: {e}")
        return pl.DataFrame()       

def run_ttest_analysis(real_df: pl.DataFrame, control_df: pl.DataFrame, n_real_total: int, n_control_total: int, p_cutoff: float) -> pl.DataFrame:
    """
    Performs group-level enrichment using a t-test on log(Odds-Ratios) with
    a robust, optimized Polars approach.
    """
    print(f"Running group comparison using T-test with a p-value cutoff of {p_cutoff}.")

    combined_df = pl.concat([
        real_df.with_columns(pl.lit("real").alias("group")),
        control_df.with_columns(pl.lit("control").alias("group"))
    ])

    analysis_df = combined_df.group_by("Trait").agg(
        pl.col("Odds-Ratio").filter(pl.col("group") == "real").alias("real_ors"),
        pl.col("Odds-Ratio").filter(pl.col("group") == "control").alias("control_ors")
    )

    def perform_ttest_on_lists(row):
        real_ors = row["real_ors"] if row["real_ors"] is not None else []
        control_ors = row["control_ors"] if row["control_ors"] is not None else []
        
        # Impute OR=1.0 for samples where the trait was not found/enriched.
        real_ors.extend([1.0] * (n_real_total - len(real_ors)))
        control_ors.extend([1.0] * (n_control_total - len(control_ors)))
        
        # If there's not enough data for a t-test, return a dictionary of NaNs
        if len(real_ors) < 2 or len(control_ors) < 2:
            return {"Statistic": np.nan, "P-value": np.nan, "Value_Real": np.nan, "Value_Control": np.nan}

        # Perform the t-test on the log-transformed Odds-Ratios.
        log_real_ors = np.log(real_ors)
        log_control_ors = np.log(control_ors)
        stat, p_value = ttest_ind(log_real_ors, log_control_ors, equal_var=False, nan_policy='omit')
        
        return {
            "Statistic": stat, 
            "P-value": p_value, 
            "Value_Real": np.mean(log_real_ors), 
            "Value_Control": np.mean(log_control_ors)
        }

    results_df = analysis_df.with_columns(
        pl.struct(["real_ors", "control_ors"]).map_elements(
            perform_ttest_on_lists,
            return_dtype=pl.Struct([
                pl.Field("Statistic", pl.Float64), 
                pl.Field("P-value", pl.Float64),
                pl.Field("Value_Real", pl.Float64),
                pl.Field("Value_Control", pl.Float64)
            ])
        ).alias("ttest_result")
    ).unnest("ttest_result").drop(["real_ors", "control_ors"]).drop_nulls()

    final_df = results_df.with_columns(
        pl.lit("T-statistic").alias("Statistic_Name"),
        pl.lit(n_real_total).alias("N_Real"),
        pl.lit(n_control_total).alias("N_Control")
    ).filter(pl.col("P-value") < p_cutoff)

    return final_df


def main():
    parser = argparse.ArgumentParser(
        description="Perform group-level enrichment analysis on CpG traits.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("--analysis-dir", type=pathlib.Path, required=True, help="Path to the module's analysis directory containing 'real_results' and 'control_results' subfolders.")
    parser.add_argument(
        "--method", 
        type=str, 
        choices=['ttest'], 
        default='ttest',
        help="""The statistical method to use for group comparison:
        'ttest': (Default) Compares the distribution of log(Odds-Ratios) using a T-test with imputation."""
    )
    parser.add_argument("--p-value-cutoff", type=float, default=0.05, help="P-value cutoff to define 'enrichment' in a single sample. Used only with --method fisher.")
    args = parser.parse_args()

    real_dir = args.analysis_dir / "real_results"
    control_dir = args.analysis_dir / "control_results"
    output_dir = args.analysis_dir

    # Count total number of sample files for imputation/frequency counts
    n_real_files = len(list(real_dir.glob('*.tsv')))
    n_control_files = len(list(control_dir.glob('*.tsv')))
    
    if n_real_files == 0:
        print(f"Error: No result files found in the real samples directory: {real_dir}")
        return
    if n_control_files == 0:
        print(f"Error: No result files found in the control samples directory: {control_dir}")
        return
        
    print(f"Found {n_real_files} real sample files and {n_control_files} control sample files.")

    # Aggregate results from all files
    real_results_df = _aggregate_results(real_dir)
    control_results_df = _aggregate_results(control_dir)

    if real_results_df.is_empty():
        print(f"Warning: No valid data found for real samples. Cannot perform T-test.")
        return
    final_df = run_ttest_analysis(real_results_df, control_results_df, n_real_files, n_control_files, args.p_value_cutoff)

    if final_df.is_empty():
        print("Warning: No traits passed the significance cutoff.")
        return

    # Save Final Results
    output_path = output_dir / f"group_comparison_dumbbell_{args.method}.tsv"
    final_df.write_csv(output_path, separator='\t')
    print(f"Successfully wrote group comparison results to {output_path}")

if __name__ == "__main__":
    main()