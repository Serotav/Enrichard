import polars as pl
import argparse
import pathlib
import re
from sys import stderr
import pyranges as pr

PROBE_ID_COL = "Probe_ID"

def annotate_cell_type(probes_df: pl.DataFrame, bed_file: pathlib.Path) -> pl.DataFrame:
    """
    Annotates probes with chromatin states for a single cell type.
    """
    eid_match = re.search(r'(E\d{3})', bed_file.name)
    if not eid_match:
        return None
    eid = eid_match.group(1)

    # Load chromatin data for the cell type
    chrom_df = pl.read_csv(
        bed_file, separator='\t', has_header=False,
        new_columns=['Chromosome', 'Start', 'End', 'State']
    )

    # Join probes with chromatin states using pyranges
    probes_pr = pr.PyRanges(probes_df.to_pandas())
    chrom_pr = pr.PyRanges(chrom_df.to_pandas())
    joined_pr = probes_pr.join(chrom_pr)

    if joined_pr.df.empty:
        return None

    # Convert back to polars and pivot
    annotated_df = pl.from_pandas(joined_pr.df)
    
    # Create a column for each state
    final_df = annotated_df.with_columns(
        pl.lit(1).cast(pl.UInt8).alias("value")
    ).pivot(
        index=PROBE_ID_COL,
        on="State",
        values="value",
        aggregate_function="max"
    ).fill_null(0)

    # Rename columns to be EID specific
    new_column_names = {col: f"{eid}_{col.split('_', 1)[1]}" for col in final_df.columns if col != PROBE_ID_COL}
    
    return final_df.rename(new_column_names)


def main():
    parser = argparse.ArgumentParser(description="Annotate probes with chromatin state data.")
    parser.add_argument('input_file', help='Probe manifest file (e.g., HM450.hg38.manifest.tsv)')
    parser.add_argument('output_file', help='Output annotated parquet file')
    parser.add_argument('cached_file', help='Output cached summary parquet file (unused, for compatibility)')
    parser.add_argument('chromatine_dir', help='Directory with chromatin state .bed files')
    parser.add_argument('metadata_file', help='EID metadata file (unused, for compatibility)')
    args = parser.parse_args()

    print("Loading probe data...", file=stderr)
    probes_df = pl.read_csv(
        args.input_file, separator='\t', null_values="NA"
    ).select(
        pl.col('Probe_ID').alias(PROBE_ID_COL),
        pl.col('CpG_chrm').alias('Chromosome'),
        pl.col('CpG_beg').alias('Start'),
        pl.col('CpG_end').alias('End')
    ).drop_nulls(['Start', 'End']).with_columns([
        pl.col('Start').cast(pl.Int64), pl.col('End').cast(pl.Int64)
    ])

    chromatine_dir = pathlib.Path(args.chromatine_dir)
    bed_files = sorted(list(chromatine_dir.glob("E*.bed.gz")))[:8] # Limiting to first 8 files

    if not bed_files:
        raise FileNotFoundError(f"No .bed.gz files starting with 'E' found in {chromatine_dir}")

    # This will be the main dataframe to join all results
    final_annotated_df = probes_df.select(PROBE_ID_COL)

    print("Annotating probes with chromatin states for each cell type...", file=stderr)
    for bed_file in bed_files:
        cell_type_df = annotate_cell_type(probes_df, bed_file)
        if cell_type_df is not None:
            final_annotated_df = final_annotated_df.join(cell_type_df, on=PROBE_ID_COL, how='left')

    # Fill any NaNs that may have resulted from the left joins
    final_annotated_df = final_annotated_df.fill_null(0)

    final_annotated_df.write_parquet(args.output_file, compression='lz4')
    print(f"Annotation complete. Saved to {args.output_file}", file=stderr)
    
    # Create a dummy cache file for compatibility with the setup script
    print("Creating dummy cache file...", file=stderr)
    pathlib.Path(args.cached_file).touch()
    print(f"Dummy cache created at {args.cached_file}", file=stderr)


if __name__ == "__main__":
    main()
