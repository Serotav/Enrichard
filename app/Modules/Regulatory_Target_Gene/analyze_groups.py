#!/usr/bin/env python3

import argparse
import pathlib
import polars as pl
import numpy as np
from scipy.stats import ttest_ind, fisher_exact

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

def run_fisher_analysis(real_df: pl.DataFrame, control_df: pl.DataFrame, n_real_total: int, n_control_total: int, p_cutoff: float) -> list:
    """
    Performs group-level enrichment by comparing the frequency of significant traits,
    filtering results *after* the Fisher's Exact Test based on the resulting p-value.
    """
    print(f"Running group comparison using Fisher's Exact Test, filtering results with a p-value cutoff of {p_cutoff}.")
    all_traits = sorted(list(set(real_df['Trait'].to_list()) | set(control_df['Trait'].to_list())))
    analysis_results = []

    for trait in all_traits:
        n_real_enriched = real_df.filter(pl.col('Trait') == trait).height
        n_control_enriched = control_df.filter(pl.col('Trait') == trait).height

        table = [
            [n_real_enriched, n_real_total - n_real_enriched],
            [n_control_enriched, n_control_total - n_control_enriched]
        ]

        if n_real_enriched == 0 and n_control_enriched == 0:
            continue

        group_or, p_value = fisher_exact(table)

        analysis_results.append({
            "Trait": trait,
            "P-value": p_value,
            "Statistic": group_or,
            "Statistic_Name": "Group_Odds_Ratio",
            "Value_Real": n_real_enriched / n_real_total if n_real_total > 0 else 0,
            "Value_Control": n_control_enriched / n_control_total if n_control_total > 0 else 0,
            "N_Real": n_real_total,
            "N_Control": n_control_total
        })
        
    filtered_results = [result for result in analysis_results if result['P-value'] < p_cutoff]
    return filtered_results

def run_ttest_analysis(real_df: pl.DataFrame, control_df: pl.DataFrame, n_real_total: int, n_control_total: int, p_cutoff: float) -> list:
    """
    Performs group-level enrichment using a t-test on log(Odds-Ratios) with imputation,
    filtering results *after* the t-test based on the resulting p-value.
    """
    print(f"Running group comparison using T-test with imputation, filtering results with a p-value cutoff of {p_cutoff}.")
    all_traits = sorted(list(set(real_df['Trait'].to_list()) | set(control_df['Trait'].to_list())))
    analysis_results = []

    for trait in all_traits:
        real_ors = real_df.filter(pl.col('Trait') == trait)['Odds-Ratio'].to_list()
        control_ors = control_df.filter(pl.col('Trait') == trait)['Odds-Ratio'].to_list()

        # Impute OR=1.0 for samples where the trait was not found/enriched
        real_ors.extend([1.0] * (n_real_total - len(real_ors)))
        control_ors.extend([1.0] * (n_control_total - len(control_ors)))

        if len(real_ors) < 2 or len(control_ors) < 2:
            continue

        # Transform data (log(1.0) = 0)
        log_real_ors = np.log(real_ors)
        log_control_ors = np.log(control_ors)

        stat, p_value = ttest_ind(log_real_ors, log_control_ors, equal_var=False, nan_policy='omit')
        
        analysis_results.append({
            "Trait": trait,
            "P-value": p_value,
            "Statistic": stat,
            "Statistic_Name": "T-statistic",
            "Value_Real": log_real_ors.mean(),
            "Value_Control": log_control_ors.mean(),
            "N_Real": len(log_real_ors),
            "N_Control": len(log_control_ors)
        })

    filtered_results = [result for result in analysis_results if result['P-value'] < p_cutoff]
    return filtered_results


def main():
    parser = argparse.ArgumentParser(
        description="Perform group-level enrichment analysis on CpG traits.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("--analysis-dir", type=pathlib.Path, required=True, help="Path to the module's analysis directory containing 'real_results' and 'control_results' subfolders.")
    parser.add_argument(
        "--method", 
        type=str, 
        choices=['fisher', 'ttest'], 
        default='fisher',
        help="""The statistical method to use for group comparison:
'fisher': Compares the frequency of enriched traits between groups using Fisher's Exact Test.
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

    analysis_results = []
    if args.method == 'fisher':
        # The fisher method can run even if one of the dataframes is empty
        analysis_results = run_fisher_analysis(real_results_df, control_results_df, n_real_files, n_control_files, args.p_value_cutoff)
    elif args.method == 'ttest':
        # The t-test method implicitly requires both dataframes to have some data to draw traits from
        if real_results_df.is_empty():
            print(f"Warning: No valid result files found in the real samples directory: {real_dir}. Cannot perform T-test.")
            return
        analysis_results = run_ttest_analysis(real_results_df, control_results_df, n_real_files, n_control_files,args.p_value_cutoff)

    if not analysis_results:
        print("Warning: No traits were eligible for statistical comparison. This can happen if no traits are enriched in any sample.")
        return

    # Save Final Results
    final_df = pl.DataFrame(analysis_results).sort("P-value")
    output_path = output_dir / f"group_comparison_dumbbell_{args.method}.tsv"
    final_df.write_csv(output_path, separator='\t')
    print(f"Successfully wrote group comparison results to {output_path}")

if __name__ == "__main__":
    main()