import pandas as pd
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
MERGE_PATH = BACKGROUND_DIR / "merge.tsv" 
PROBE_ID_COL = "Probe_ID"

USER_CUSTOM_BACKGROUND_NAME = os.getenv("USER_CUSTOM_BACKGROUND_NAME")
USER_SAMPLE_NAME = os.getenv("USER_SAMPLE_NAME")
USER_MODULE_OUTPUT = os.getenv("USER_MODULE_OUTPUT")
OUTPUT_FILE_NAME = str(APP_DIR.name) + ".tsv"
LOOK_UP_TABLE = APP_DIR / "lookup.tsv"

def merge_background():
    
    source_files = [
        file_path for file_path in BACKGROUND_DIR.glob("*.tsv")
        if file_path.name != MERGE_PATH.name # Exclude the merge file itself, WE SHOULD NOT BE HERE IF IT EXISTS
    ]

    if not source_files:
        print(f"No TSV files to merge in {BACKGROUND_DIR} directory.", file=sys.stderr)
        return pl.DataFrame() # Return empty DataFrame if no files

    # ALWAYS READ WITH THE RIGHT TYPES OR INT64 INSTREAD OF UINT64 WILL CRASH US  
    first_file_path = source_files[0]
    with open(first_file_path, 'r') as f:
        header = f.readline().strip().split('\t')
    
    schema_overrides = {col: pl.UInt64 for col in header if col.startswith("bitset_")}

    lazy_frames = [
        pl.scan_csv(file_path, separator='\t', schema_overrides=schema_overrides)
        for file_path in source_files
    ]

    #  Concatenate all lazy frames, drop duplicates, and then collect the result.
    merged_df = pl.concat(lazy_frames).unique(subset=[PROBE_ID_COL], keep='first').collect()

    merged_df.write_csv(MERGE_PATH, separator='\t')
    return merged_df

def custom_background(user_background_file: Path) -> pl.DataFrame:
    if not MERGE_PATH.exists():
        merge_df = merge_background()
    else: 
        with open(MERGE_PATH, 'r') as f:
            header = f.readline().strip().split('\t')
        
        # ALWAYS READ WITH THE RIGHT TYPES OR INT64 INSTREAD OF UINT64 WILL CRASH US  
        schema_overrides = { col: pl.UInt64 for col in header if col.startswith("bitset_")}
        merge_df = pl.read_csv(MERGE_PATH, separator='\t', schema_overrides=schema_overrides)
    
    user_probes_df = pl.read_csv(user_background_file, has_header=False).rename({"column_1": PROBE_ID_COL})
    custom_annotated_df = merge_df.join(user_probes_df, on=PROBE_ID_COL, how="inner")
    
    return custom_annotated_df


def perform_enrichment(background_bitset_df: pl.DataFrame, sample_probes_df: pl.DataFrame) -> pl.DataFrame:
    """
    Performs Fisher's exact test using a Polars backend for fast counting.
    """
    #  Preparation 
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

    # Iterate and Test Each Trait
    # TODO: rewrite this to use Polars expressions instead of iterating over rows
    for row in lookup_df.iter_rows(named=True):
        trait_name, bit_pos = row['trait'], row['index']
        
        # Determine the correct bitset column and the bitmask to check for this trait
        chunk_id = bit_pos // BITS_PER_CHUNK
        bitmask_to_check = 1 << (bit_pos % BITS_PER_CHUNK)
        col_to_check = f"bitset_{chunk_id}"

        # This expression checks if the trait's bit is set.
        has_trait_expr = (pl.col(col_to_check) & bitmask_to_check) > 0

        # Build the 2x2 contingency table using Polars expressions
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
        
        # Run Fisher's Exact Test 
        odds_ratio, p_value = fisher_exact([[a, b], [c, d]], alternative='greater')
        
        # Store significant results 
        results_list.append({
            "Trait": trait_name,
            "P-Value": p_value,
            "Odds-Ratio": odds_ratio,
            "a": a, "b": total_sample_size - a, "c": c, "d": total_background_only_size - c,
        })
            
    # Create and return the final DataFrame 
    if not results_list:
        return pl.DataFrame()
        
    return pl.DataFrame(results_list).sort("P-Value")

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
        print("Applying no multiple testing correction.", file=sys.stderr)
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
    
    print(f"Applied '{method}' correction. Found {len(significant_results)} significant results.", file=sys.stderr)
    
    return significant_results.sort("P-adj")

def get_background_df(background_name: str, user_dir: Path):
    if background_name == USER_CUSTOM_BACKGROUND_NAME:
        print("Processing custom background...", file=sys.stderr)
        # The user's custom background file was saved in their session directory
        user_background_file = user_dir / USER_CUSTOM_BACKGROUND_NAME
        return custom_background(user_background_file)
    
    # Find a file in BACKGROUND_DIR that starts with args.background_name
    background_dir = pathlib.Path(BACKGROUND_DIR)
    matching_files = list(background_dir.glob(f"{background_name}*"))
    if not matching_files or len(matching_files) > 1:
        print(f"Error: No file found in {background_dir} starting with '{background_name}' or multiple found files found:[{len(matching_files)}]", file=sys.stderr)
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
    parser.add_argument("--correction", default='none',  help="Multiple testing correction method to apply.")

    args = parser.parse_args()
    user_dir = Path(args.user_dir)

    # Load Data 
    try:
        background_df = get_background_df(args.background_name, user_dir)
        sample_df = pl.read_csv(Path(args.user_dir) / USER_SAMPLE_NAME, has_header=False).rename({ "column_1": PROBE_ID_COL })
    except Exception as e:
        print(f"Error reading background or sample {e}", file=sys.stderr)
        sys.exit(1)

    # Run Analysis 
    print("Starting enrichment analysis...", file=sys.stderr)
    raw_results_df = perform_enrichment(background_df, sample_df)
    final_results_df = apply_correction(raw_results_df, args.correction, args.p_value)

    # Save Results 
    # Check if output directory exists, create if not
    output_dir = pathlib.Path(args.user_dir) / USER_MODULE_OUTPUT / APP_DIR.name
    output_dir.mkdir(parents=True, exist_ok=True)
    destination = output_dir/ OUTPUT_FILE_NAME
    
    final_results_df.write_csv(destination, separator="\t")
    print(f"Significant results saved to: {destination}")


if __name__ == "__main__":
    main()