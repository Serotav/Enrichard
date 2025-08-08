### Git lfs
This module needs git lfs to download the pre annotated data
```bash
sudo pacman -S git-lfs
git lfs install
git lfs pull
```
not on arch? *skill issue*.


this is how the brackground files were annotated:
```bash
#!/bin/bash
set -e

# This script downloads histone mark bigWig files (E001-E129).
# It first verifies that all E001 files exist on the remote server before starting.
SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
BASE_URL="https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed"
HISTONE_FOLDER="$SCRIPT_DIR/histone_data"
HISTONE_MARKS=(
    "DNase" "H2A.Z" "H3K27ac" "H3K27me3" "H3K4me1" "H3K4me2"
    "H3K4me3" "H3K79me2" "H3K9ac" "H3K9me3" "H4K20me1"
)
CELL_TYPE_COUNT=129

# Remote File Verification (using E001 as a test case) 
echo "Phase 1: Verifying remote file availability using E001..."
for mark in "${HISTONE_MARKS[@]}"; do
    filename="E001-${mark}.imputed.pval.signal.bigwig"
    url="$BASE_URL/$mark/$filename"
    echo -n "  - Checking for $mark... "
    if ! wget --spider -q "$url"; then
        echo "FAILED."
        echo "Error: Remote file not found at $url"
        echo "Aborting script to prevent partial downloads."
        exit 1
    fi
    echo "OK."
done
echo "Remote file verification successful."
echo "---"

# Local Data Check 
echo "Phase 2: Checking for existing local data..."
if [ -d "$SCRIPT_DIR/DNase" ] && [ -n "$(ls -A "$SCRIPT_DIR/DNase")" ]; then
    echo "Histone directories already exist and contain data. Skipping download."
    exit 0
else
    echo "No existing data found. Proceeding with download."
fi
echo "---"

echo "Creating subdirectories for each histone mark..."
mkdir -p $HISTONE_FOLDER
cd $HISTONE_FOLDER
for mark in "${HISTONE_MARKS[@]}"; do
    mkdir -p "$SCRIPT_DIR/$mark"
done
echo "Subdirectories created."

# Loop from 1 to 129 to download all files
for i in $(seq 1 $CELL_TYPE_COUNT); do
    # Format the number to be three digits (1 -> 001)
    eid_num=$(printf "%03d" $i)
    eid="E${eid_num}"
    
    echo "Downloading data for cell type: $eid"
    
    for mark in "${HISTONE_MARKS[@]}"; do
        filename="${eid}-${mark}.imputed.pval.signal.bigwig"
        url="$BASE_URL/$mark/$filename"
        target_dir="$SCRIPT_DIR/$mark"
        
        # Check if the specific file already exists before downloading
        if [ -f "$target_dir/$filename" ]; then
            echo "  - Found ${filename}, skipping."
        else
            echo "  - Downloading ${filename} into $mark/ ..."
            # Use wget to download quietly, but check for server errors first
            if wget --spider -q "$url"; then
                wget -q -P "$target_dir" "$url"
            else
                # This case should be rare since we verified E001, but it's good practice
                echo "    ...Warning: File not found on server at $url. Skipping."
            fi
        fi
    done
    echo "---"
done

echo "Download script finished successfully."
```

```python
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

def get_sort_key(file_path):
    path = pathlib.Path(file_path)
    enzyme = path.parent.name
    eid_match = re.search(r'(E\d{3})', path.name)

    if eid_match:
        return (eid_match.group(1), enzyme)

    # if no EID, we're cooked
    return ('Z999', enzyme)


if __name__ == "__main__":
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

    bw_files.sort(key=get_sort_key)

    print("File processing order is locked in, periodt.", file=stderr)
    for f in bw_files:
        print(f"  -> {f}", file=stderr)

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

```