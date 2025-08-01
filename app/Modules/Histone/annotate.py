

import polars as pl
import argparse
import pathlib
import re
from sys import stderr
import pyBigWig
import glob
import numpy as np

PROBE_ID_COL = "Probe_ID"
# The threshold to determine if an enzyme binds. The data is in -log10(p-value)
# format, so we use a threshold of -log10(0.05), which is ~1.3.
VALUE_THRESHOLD = -np.log10(0.05)

def main():
    parser = argparse.ArgumentParser(description="Annotate probes with histone data from BigWig files.")
    parser.add_argument('input_file', help='Probe manifest file (e.g., HM450.hg38.manifest.tsv)')
    parser.add_argument('output_file', help='Output annotated parquet file')
    parser.add_argument('cached_file', help='Output cached summary parquet file (unused, for compatibility)')
    parser.add_argument('histone_dir', help='Directory with histone BigWig files')
    args = parser.parse_args()

    print("Loading probe data...", file=stderr)
    probes_df = pl.read_csv(
        args.input_file, separator='\t', null_values="NA"
    ).select(
        pl.col('Probe_ID').alias(PROBE_ID_COL),
        pl.col('CpG_chrm').alias('Chromosome'),
        pl.col('CpG_beg').alias('Start'),
    ).drop_nulls(['Start']).with_columns([
        pl.col('Start').cast(pl.Int64)
    ])

    histone_dir = pathlib.Path(args.histone_dir)
    bw_files = glob.glob(f"{histone_dir}/**/*.bigwig", recursive=True)

    if not bw_files:
        raise FileNotFoundError(f"No .bigwig files found in {histone_dir}")

    final_annotated_df = probes_df.select(PROBE_ID_COL)

    print(f"Found {len(bw_files)} BigWig files to process.", file=stderr)

    for i, bw_file_path in enumerate(bw_files):
        bw_file = pathlib.Path(bw_file_path)
        enzyme = bw_file.parent.name
        
        eid_match = re.search(r'(E\d{3})', bw_file.name)
        if not eid_match:
            print(f"Could not extract EID from {bw_file.name}, skipping.", file=stderr)
            continue
        
        eid = eid_match.group(1)
        column_name = f"{eid}_{enzyme}"
        
        print(f"Processing file {i+1}/{len(bw_files)}: {bw_file.name} -> {column_name}", file=stderr)

        try:
            bw = pyBigWig.open(str(bw_file))
        except Exception as e:
            print(f"Error opening BigWig file {bw_file}: {e}", file=stderr)
            continue

        annotations = []
        for chrom, start in probes_df.select("Chromosome", "Start").iter_rows():
            if chrom is None or start is None:
                annotations.append(0)
                continue
            
            try:
                # Query the mean value over a 50bp window starting from the probe's start
                mean_val = bw.stats(chrom, int(start), int(start) + 50, type='mean')
                if mean_val is not None and mean_val[0] is not None and mean_val[0] > VALUE_THRESHOLD:
                    annotations.append(1)
                else:
                    annotations.append(0)
            except RuntimeError:
                # This can happen if the chromosome is not in the BigWig file
                annotations.append(0)

        bw.close()
        
        final_annotated_df = final_annotated_df.with_columns(
            pl.Series(name=column_name, values=annotations, dtype=pl.UInt8)
        )

    final_annotated_df.write_parquet(args.output_file, compression='lz4')
    print(f"Annotation complete. Saved to {args.output_file}", file=stderr)

    # Create a dummy cache file for compatibility with the setup script
    print("Creating dummy cache file...", file=stderr)
    pathlib.Path(args.cached_file).touch()
    print(f"Dummy cache created at {args.cached_file}", file=stderr)


if __name__ == "__main__":
    main()

