#!/usr/bin/env python3

import argparse
import pathlib
import polars as pl
import numpy as np
from scipy.stats import ttest_ind

def _aggregate_and_melt_heatmaps(directory: pathlib.Path) -> pl.DataFrame:
    """
    Reads all wide-format TSV heatmap files in a directory, melts them
    into a long format, and concatenates them into a single DataFrame.
    """
    all_long_dfs = []

    for file_path in directory.glob('*.tsv'):
        try:
            # The first column is the Cell_Type index
            df = pl.read_csv(file_path, separator='\t')
            if df.is_empty():
                print(f"Warning: Skipping empty file: {file_path}")
                continue

            # Melt the wide-format data into a long format
            # The first column name is used as the id_vars
            id_col_name = df.columns[0]
            long_df = df.melt(
                id_vars=[id_col_name],
                variable_name="State",
                value_name="Odds-Ratio"
            ).rename({id_col_name: "Cell_Type"})
            
            all_long_dfs.append(long_df)

        except pl.exceptions.NoDataError:
            print(f"Warning: Skipping empty or malformed file: {file_path}")
            continue
        except Exception as e:
            print(f"An unexpected error occurred while processing {file_path}: {e}")
            continue
    
    if not all_long_dfs:
        return pl.DataFrame()
        
    return pl.concat(all_long_dfs)


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
    
    # Define a specific output filename for this module's analysis
    output_filename = "group_comparison_heatmap_ttest.tsv"
    output_path = output_dir / output_filename
    
    # Count total number of sample files for imputation
    n_real_total = len(list(real_dir.glob('*.tsv')))
    n_control_total = len(list(control_dir.glob('*.tsv')))

    if n_real_total == 0 or n_control_total == 0:
        print("Error: 'real_results' and 'control_results' directories must both contain at least one file.", file=sys.stderr)
        return

    # Aggregate and transform data from both groups
    print("Aggregating results from real samples...")
    real_long_df = _aggregate_and_melt_heatmaps(real_dir)
    print("Aggregating results from control samples...")
    control_long_df = _aggregate_and_melt_heatmaps(control_dir)

    if real_long_df.is_empty() and control_long_df.is_empty():
        print("Warning: No valid data found in either real or control directories.")
        return

    # Get all unique (Cell_Type, State) pairs to iterate over
    all_pairs = pl.concat([
        real_long_df.select(["Cell_Type", "State"]),
        control_long_df.select(["Cell_Type", "State"])
    ]).unique()

    print(f"Found {len(all_pairs)} unique (Cell_Type, State) pairs to test...")
    analysis_results = []

    # --- Main Analysis Loop ---
    for row in all_pairs.iter_rows(named=True):
        cell_type, state = row['Cell_Type'], row['State']

        # Get all Odds-Ratios for the current pair from both groups
        real_ors = real_long_df.filter(
            (pl.col('Cell_Type') == cell_type) & (pl.col('State') == state)
        )['Odds-Ratio'].to_list()
        
        control_ors = control_long_df.filter(
            (pl.col('Cell_Type') == cell_type) & (pl.col('State') == state)
        )['Odds-Ratio'].to_list()

        # Impute OR=1.0 for samples where the pair was not found
        real_ors.extend([1.0] * (n_real_total - len(real_ors)))
        control_ors.extend([1.0] * (n_control_total - len(control_ors)))
        
        # Log-transform data; log(1.0) = 0
        log_real_ors = np.log(real_ors)
        log_control_ors = np.log(control_ors)

        # Perform Welch's t-test
        stat, p_value = ttest_ind(log_real_ors, log_control_ors, equal_var=False, nan_policy='omit')
        
        analysis_results.append({
            "Cell_Type": cell_type,
            "State": state,
            "P-value": p_value,
            "T-statistic": stat,
            "Mean_logOR_Real": np.mean(log_real_ors),
            "Mean_logOR_Control": np.mean(log_control_ors)
        })

    if not analysis_results:
        print("Warning: No results were generated from the analysis.")
        return

    # Convert to DataFrame and apply the final p-value filter
    results_df = pl.DataFrame(analysis_results)
    final_df = results_df.filter(pl.col("P-value") < args.p_value_cutoff).sort("P-value")

    if final_df.is_empty():
        print(f"Warning: No (Cell_Type, State) pairs passed the significance threshold of p < {args.p_value_cutoff}.")
        # Create an empty file to signify completion
        output_path.touch()
        return

    # Save the final long-format results
    final_df.write_csv(output_path, separator='\t')
    print(f"Successfully wrote group comparison results for chromatin to {output_path}")

if __name__ == "__main__":
    main()