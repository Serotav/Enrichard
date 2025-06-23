import pandas as pd
import polars as pl
from scipy.stats import fisher_exact
import sys
import argparse
import os
import glob
import pathlib
from pathlib import Path
APP_DIR = pathlib.Path(__file__).parent
BACKGROUND_DIR  = APP_DIR/"background"
MERGE = BACKGROUND_DIR/"merge.tsv"
PROBE_ID_COL = "Probe_ID"

USER_CUSTOM_BACKGROUND_NAME = os.getenv("USER_CUSTOM_BACKGROUND_NAME")
USER_SAMPLE_NAME = os.getenv("USER_SAMPLE_NAME")
USER_MODULE_OUTPUT = os.getenv("USER_MODULE_OUTPUT")
OUTPUT_FILE_NAME = str(APP_DIR.name) + ".tsv"
LOOK_UP_TABLE = APP_DIR / "lookup.tsv"

def merge_background():
    print(f'Merging TSV files in {BACKGROUND_DIR}', file=sys.stderr)
    tsv_files = glob.glob(os.path.join(BACKGROUND_DIR, "*.tsv"))

    if not tsv_files:
        print(f"No TSV files found in {BACKGROUND_DIR} directory.", file=sys.stderr)
        return

    all_data = []
    for file_path in tsv_files:
        df = pd.read_csv(file_path, sep='\t')
        all_data.append(df)

    merged_df = pd.concat(all_data, ignore_index=True)
    print(f"Total rows before deduplication: {len(merged_df)}", file=sys.stderr)
    
    # It's better to drop duplicates based on the Probe_ID to ensure one annotation per probe
    merged_df = merged_df.drop_duplicates(subset=['Probe_ID'])
    print(f"Total rows after deduplication on Probe_ID: {len(merged_df)}", file=sys.stderr)

    save_path = pathlib.Path(BACKGROUND_DIR) / MERGE
    merged_df.to_csv(save_path, sep='\t', index=False)
    print(f"Merged file saved to: {save_path}", file=sys.stderr)
    return merged_df

def custom_background(background_file):

    if not pathlib.Path(BACKGROUND_DIR/MERGE).exists():
        merge = merge_background()
    else: 
        print(f"File {os.path.join(output_dir, file_name)} already exists, skipping merge.", file=sys.stderr)
        merge = pd.read_csv(os.path.join(output_dir, file_name), sep='\t')
    
    # The user's uploaded background is a simple list of probes, likely with no header.
    user_background = pd.read_csv(background_file, header=None)
    user_background.columns = ['Probe_ID']

    # Filter the large merged annotation file to only include the user's probes.
    filtered_merge = merge[merge['Probe_ID'].isin(user_background['Probe_ID'])]
    
    # Overwrite the original uploaded file with this new, filtered, and annotated version.
    filtered_merge.to_csv(background_file, sep='\t', index=False)
    print(f"Created custom annotated background file at: {background_file}", file=sys.stderr)


def perform_enrichment(background_bitset_df: pl.DataFrame, sample_probes_df: pl.DataFrame, p_value_threshold: float) -> pl.DataFrame:
    """
    Performs Fisher's exact test using a Polars backend for fast counting.
    """
    # --- Step 1: Preparation ---
    try:
        lookup_df = pl.read_csv(LOOK_UP_TABLE, separator='\t')
    except FileNotFoundError:
        print(f"Error: Lookup table not found at {LOOK_UP_TABLE}", file=sys.stderr)
        sys.exit(1)

    # Pre-calculate which probes in the background are part of the user's sample.
    background_df = background_bitset_df.with_columns(
        is_in_sample=pl.col(PROBE_ID_COL).is_in(sample_probes_df[PROBE_ID_COL].unique())
    )
    
    # Pre-calculate total counts for efficiency inside the loop
    total_sample_size = background_df.filter(pl.col("is_in_sample")).height
    total_background_only_size = background_df.height - total_sample_size
    
    BITS_PER_CHUNK = 64
    results_list = []

    # --- Step 2: Iterate and Test Each Trait ---
    for row in lookup_df.iter_rows(named=True):
        trait_name, bit_pos = row['trait'], row['index']
        
        # Determine the correct bitset column and the bitmask to check for this trait
        chunk_id = bit_pos // BITS_PER_CHUNK
        bitmask_to_check = 1 << (bit_pos % BITS_PER_CHUNK)
        col_to_check = f"bitset_{chunk_id}"

        # FIX 2: Use the `&` operator for bitwise AND instead of the .bitwise_and() method.
        # This expression checks if the trait's bit is set.
        has_trait_expr = (pl.col(col_to_check) & bitmask_to_check) > 0

        # Efficiently build the 2x2 contingency table using Polars expressions
        counts = background_df.select(
            # 'a': Probes IN sample that HAVE the trait
            a=has_trait_expr.filter(pl.col("is_in_sample")).sum(),
            # 'c': Probes NOT in sample that HAVE the trait
            c=has_trait_expr.filter(~pl.col("is_in_sample")).sum()
        ).row(0)
        
        a, c = counts[0], counts[1]

        # If no probes have this trait, skip the test
        if a + c == 0:
            continue
        
        # Calculate 'b' and 'd' from the totals
        b = total_sample_size - a
        d = total_background_only_size - c
        
        # --- Run Fisher's Exact Test ---
        odds_ratio, p_value = fisher_exact([[a, b], [c, d]], alternative='greater')
        
        # --- Store significant results ---
        if p_value < p_value_threshold:
            results_list.append({
                "Trait": trait_name,
                "P-Value": p_value,
                "Odds-Ratio": odds_ratio,
                "a": a, "b": b, "c": c, "d": d,
            })
            
    # --- Step 3: Create and return the final DataFrame ---
    if not results_list:
        return pl.DataFrame()
        
    return pl.DataFrame(results_list).sort("P-Value")

def get_brackground_df(background:str):
    if background == USER_CUSTOM_BACKGROUND_NAME:
        print("Detected custom background path. Not implemented", file=sys.stderr)
        exit(1)
        return custom_background(args.background_name)

    # Find a file in BACKGROUND_DIR that starts with args.background_name
    background_dir = pathlib.Path(BACKGROUND_DIR)
    matching_files = list(background_dir.glob(f"{background}*"))
    if not matching_files or len(matching_files) > 1:
        print(f"Error: No file found in {background_dir} starting with '{background}' or multiple found files found:[{len(matching_files)}]", file=sys.stderr)
        sys.exit(1)
    
    # Read with the right types
    file_path = matching_files[0]
    with open(file_path, 'r') as f:
        header = f.readline().strip().split('\t')
    schema_overrides = {col: pl.UInt64 for col in header[1:]} 
    df = pl.read_csv(
        file_path,
        separator='\t',
        null_values="NA",
        schema_overrides=schema_overrides
    )
    return df


def main():

    parser = argparse.ArgumentParser(description="Perform Fisher's exact test for enrichment analysis on CpG sites.")
    parser.add_argument("--user_dir", help="Path to the user directory.")
    parser.add_argument("--background_name", help="Name of the background annotation file.")
    parser.add_argument("--p_value", type=float, default=0.05, help="P-value threshold for significance.")
    args = parser.parse_args()

    # Load Data 
    try:
        background_df = get_brackground_df(args.background_name)
        sample_df = pl.read_csv(Path(args.user_dir) / USER_SAMPLE_NAME, has_header=False).rename({ "column_1": PROBE_ID_COL })
    except Exception as e:
        print(f"Error reading background or sample {e}", file=sys.stderr)
        sys.exit(1)

    # --- Run Analysis ---
    print("Starting enrichment analysis...", file=sys.stderr)
    
    results_df = perform_enrichment(background_df, sample_df, args.p_value)

    # --- Save Results ---
    # Check if output directory exists, create if not
    output_dir = pathlib.Path(args.user_dir) / USER_MODULE_OUTPUT / APP_DIR.name
    output_dir.mkdir(parents=True, exist_ok=True)
    destination = output_dir/ OUTPUT_FILE_NAME
    
    results_df = results_df.sort("P-Value")
    results_df.write_csv(destination, separator="\t")
    print(f"Significant results saved to: {destination}")


if __name__ == "__main__":
    main()