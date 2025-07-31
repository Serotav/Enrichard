#!/bin/bash
set -e

# This script downloads histone mark bigWig files (E001-E129).
# It first verifies that all E001 files exist on the remote server before starting.

# --- Configuration ---
SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
BASE_URL="https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed"
HISTONE_FOLDER="$SCRIPT_DIR/histone_data"
HISTONE_MARKS=(
    "DNase" "H2A.Z" "H3K27ac" "H3K27me3" "H3K4me1" "H3K4me2"
    "H3K4me3" "H3K79me2" "H3K9ac" "H3K9me3" "H4K20me1"
)
CELL_TYPE_COUNT=129

# --- Phase 1: Remote File Verification (using E001 as a test case) ---
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

# --- Phase 2: Local Data Check ---
echo "Phase 2: Checking for existing local data..."
# Check if the first histone directory exists and has files
if [ -d "$SCRIPT_DIR/DNase" ] && [ -n "$(ls -A "$SCRIPT_DIR/DNase")" ]; then
    echo "Histone directories already exist and contain data. Skipping download."
    exit 0
else
    echo "No existing data found. Proceeding with download."
fi
echo "---"

# --- Phase 3: Download and Organize ---
echo "Phase 3: Starting download and organization..."

# Create subdirectories for each histone mark
echo "Creating subdirectories for each histone mark..."
mkdir -p $HISTONE_FOLDER
cd $HISTONE_FOLDER
for mark in "${HISTONE_MARKS[@]}"; do
    mkdir -p "$SCRIPT_DIR/$mark"
done
echo "Subdirectories created."

# Loop from 1 to 129 to download all files
for i in $(seq 1 $CELL_TYPE_COUNT); do
    # Format the number to be three digits (e.g., 1 -> 001, 12 -> 012, 129 -> 129)
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
