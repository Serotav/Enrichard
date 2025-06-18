import pandas as pd
from scipy.stats import fisher_exact
import numpy as np
from collections import defaultdict
import sys
import argparse
import os
import glob

def merge_background(background_dir, output_dir,file_name):
    print(f'Merging TSV files in {background_dir}', file=sys.stderr)
    tsv_files = glob.glob(os.path.join(background_dir, "*.tsv"))

    if not tsv_files:
        print(f"No TSV files found in {background_dir} directory. Current working directory: {os.getcwd()}", file=sys.stderr)
        return

    # Read and concatenate all TSV files
    all_data = []
    for file_path in tsv_files:
        df = pd.read_csv(file_path, sep='\t')
        all_data.append(df)

    # Concatenate all dataframes
    merged_df = pd.concat(all_data, ignore_index=True)
    print(f"Total rows before deduplication: {len(merged_df)}", file=sys.stderr)

    # Remove duplicates
    merged_df = merged_df.drop_duplicates()
    print(f"Total rows after deduplication: {len(merged_df)}", file=sys.stderr)

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Save merged file
    output_path = os.path.join(output_dir, file_name)
    merged_df.to_csv(output_path, sep='\t', index=False)
    print(f"Merged file saved to: {output_path}", file=sys.stderr)
    return merged_df

def custom_background(background_file):
    background_dir = "background"
    output_dir = "background/merge"
    file_name = "merge.tsv"
    
    if not os.path.exists(os.path.join(output_dir, file_name)):
        merge = merge_background(background_dir, output_dir,file_name)
    else: 
        print(f"File {os.path.join(output_dir, file_name)} already exists, skipping merge.", file=sys.stderr)
        merge = pd.read_csv(os.path.join(output_dir, file_name), sep='\t')
    
    user_background = pd.read_csv(background_file, sep='\t')

    # Rename the first column of user_background to 'Probe_ID'
    user_background.columns = ['Probe_ID'] + list(user_background.columns[1:])
    # Filter merge to only include rows where Probe_ID is in user_background
    filtered_merge = merge[merge['Probe_ID'].isin(user_background['Probe_ID'])]
    filtered_merge.to_csv(background_file, sep='\t', index=False)

def main():
    print(sys.argv, file=sys.stderr)
    # Here we set up the sample and background files, the output file, and the cols that containts the traits
    parser = argparse.ArgumentParser(description="Perform Fisher's exact test for enrichment analysis on CpG sites.")
    parser.add_argument("background_file", help="Path to the background annotation file (e.g., MSA_parsed.csv).")
    parser.add_argument("sample_file", help="Path to the sample file containing probe IDs (one per line, no header).")
    parser.add_argument("output_file", help="Path to save the significant enrichment results (CSV format).")
    parser.add_argument("--p_value_threshold", type=float, default=0.05, help="P-value threshold for significance (default: 0.05).")
    parser.add_argument("--cols_contain", type=str, default='hs', help="Substring to identify trait columns (default: 'hs').")


    COLNAME = "IlmnID" # "probeID"
    args = parser.parse_args()
    if '/sample/' in args.background_file:
        custom_background(args.background_file)

    print(f"ok")
    exit(0)

    # Check if the files exist and read them
    # Also, we will merge the sample with the background on probeID, aka we are dropping the probes that are not in the background
    print(f"Reading background file: {args.background_file}, sample {args.sample_file}", file=sys.stderr)
    try:
        background = pd.read_csv(args.background_file, sep='\t')
        background = background[background["ENTREZID"].notna()]
        sample = pd.read_csv(args.sample_file, header=None)
        sample.columns = [COLNAME] # Assume single column file with probe IDs
        sample = pd.merge(sample, background, on=COLNAME, how="inner")
    except Exception as e:
        print(f"Error reading or merging file", file=sys.stderr)
        sys.exit(1)

    # We get the trait columns from the background file based on the substring provided
    trait_cols = [col for col in background.columns if args.cols_contain in col.lower()] 
    results = []

    print('args',args.cols_contain,trait_cols, file=sys.stderr)
    

    # !!!!!!!! :1 is for testing !!!!!!!!!!!!!!!!
    for trait_col in trait_cols[:2]:
        # Prepare background data for the current trait column
        back_trait_map = defaultdict(lambda: 0)
        multiple_back_entries = []
        valid_back_entries = background[trait_col].dropna()
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
        valid_sample_entries = sample[trait_col].dropna()
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
                             "Trait Column": trait_col,
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
                      print(f"Skipping Fisher test for trait '{trait}' in column '{trait_col}' due to error: {e}", file=sys.stderr)
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