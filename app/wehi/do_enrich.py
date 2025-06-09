import pandas as pd
from scipy.stats import fisher_exact
import numpy as np
from collections import defaultdict
import sys
import argparse

def count_in_multiple(multiple:list[str], target:str):
    count = 0
    for entry in multiple:
        if target in entry:
            count += 1
    return count

def main():

    parser = argparse.ArgumentParser(description="Perform Fisher's exact test for enrichment analysis on CpG sites.")
    parser.add_argument("background_file", help="Path to the background annotation file (e.g., MSA_parsed.csv).")
    parser.add_argument("sample_file", help="Path to the sample file containing probe IDs (one per line, no header).")
    parser.add_argument("output_file", help="Path to save the significant enrichment results (CSV format).")
    parser.add_argument("--p_value_threshold", type=float, default=0.05, help="P-value threshold for significance (default: 0.05).")
    parser.add_argument("--cols_contain", type=str, default='human', help="Substring to identify trait columns (default: 'human').")

    args = parser.parse_args()

    try:
        background = pd.read_csv(args.background_file, sep='\t')
        background = background[background["ENTREZID"].notna()]
        sample = pd.read_csv(args.sample_file, header=None)
        sample.columns = ["probeID"] # Assume single column file with probe IDs
        sample = pd.merge(sample, background, on="probeID", how="inner")
    except Exception as e:
        print(f"Error reading or merging file", file=sys.stderr)
        sys.exit(1)

   
    trait_cols = [col for col in background.columns if args.cols_contain in col.lower()] 
    results = []


    # !!!!!!!! :1 is for testing !!!!!!!!!!!!!!!!
    for col in trait_cols[:1]:
        # Prepare background data for the current trait column
        back_trait_map = defaultdict(lambda: 0)
        multiple_back_entries = []
        valid_back_entries = background[col].dropna()
        for entry in valid_back_entries:
            if ';' in entry:
                multiple_back_entries.extend(entry.split(';')) # Store individual traits
            else:
                back_trait_map[entry] += 1
        
        # Count unique traits in background, handling single and multiple entries
        unique_back_traits = set(back_trait_map.keys()) | set(multiple_back_entries)
        background_total_annotations = sum(back_trait_map.values()) + len(multiple_back_entries)


        # Prepare sample data for the current trait column
        sample_trait_map = defaultdict(lambda: 0)
        multiple_sample_entries = []
        valid_sample_entries = sample[col].dropna()
        for entry in valid_sample_entries:
             if ';' in entry:
                 multiple_sample_entries.extend(entry.split(';')) # Store individual traits
             else:
                 sample_trait_map[entry] += 1
        
        sample_total_annotations = sum(sample_trait_map.values()) + len(multiple_sample_entries)


        # Perform Fisher's exact test for each unique trait found in the background
        for trait in unique_back_traits:
             # Counts for contingency table
             sample_has_trait = sample_trait_map[trait] + multiple_sample_entries.count(trait)
             sample_no_trait = sample_total_annotations - sample_has_trait
             
             back_has_trait = back_trait_map[trait] + multiple_back_entries.count(trait)
             back_no_trait = background_total_annotations - back_has_trait

             # Ensure counts are non-negative (can happen with edge cases/data issues)
             sample_has_trait = max(0, sample_has_trait)
             sample_no_trait = max(0, sample_no_trait)
             back_has_trait = max(0, back_has_trait)
             back_no_trait = max(0, back_no_trait)


             # Create contingency table
             # [[sample_with, background_with], [sample_without, background_without]]
             contingency_table = np.array([[sample_has_trait, back_has_trait],
                                         [sample_no_trait, back_no_trait]])

             # Perform the test only if there's data to test
             if np.sum(contingency_table) > 0 and sample_has_trait + back_has_trait > 0 : # Check if trait present in sample or background
                 try:
                     odds_ratio, p_value = fisher_exact(contingency_table, alternative='greater') # Test for enrichment
                     
                     # Check for significance
                     if p_value < args.p_value_threshold:
                         results.append({
                             "Trait Column": col,
                             "Trait": trait,
                             "Odds Ratio": odds_ratio,
                             "P-Value": p_value,
                             "Sample with Trait": sample_has_trait,
                             "Sample Total Annotations": sample_total_annotations,
                             "Background with Trait": back_has_trait,
                             "Background Total Annotations": background_total_annotations
                         })
                 except ValueError as e:
                      # This can happen if the table contains non-integer/negative values (shouldn't with the max(0,...) checks)
                      # or if the sum of rows/columns is zero in some edge cases.
                      print(f"Skipping Fisher test for trait '{trait}' in column '{col}' due to error: {e}", file=sys.stderr)
                      print(f"Contingency Table: {contingency_table.tolist()}", file=sys.stderr)
    if results:
        results_df = pd.DataFrame(results)
        results_df = results_df.sort_values(by="P-Value")
        try:
            results_df.to_csv(args.output_file, index=False)
            print(f"Significant results saved to: {args.output_file}")
        except Exception as e:
             print(f"Error saving results to {args.output_file}: {e}", file=sys.stderr)
             sys.exit(1)
    else:

        print(f"No significant enrichment found (p < {args.p_value_threshold}). No output file generated.")
        # Create an empty file to signal completion if needed by the app
        try:
            open(args.output_file, 'w').close()
            print(f"Empty results file created at: {args.output_file}")
        except Exception as e:
            print(f"Error creating empty results file {args.output_file}: {e}", file=sys.stderr)



if __name__ == "__main__":
    main()