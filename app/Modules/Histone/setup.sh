#!/bin/bash
set -e
SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"

# --- Configuration ---
PRECOMPUTED_URL_FILE="$SCRIPT_DIR/precomputed_url.txt"
RAW_DIR="$SCRIPT_DIR/raw_histone"
BACKGROUND_DIR="$SCRIPT_DIR/background"
CACHED_DIR="$SCRIPT_DIR/cached"
ANNOTATE_SCRIPT_PY="$SCRIPT_DIR/annotate.py"

# --- Create necessary directories ---
mkdir -p "$BACKGROUND_DIR" "$CACHED_DIR"

# --- Helper Functions ---

download_and_organize_raw_data() {
    echo "--- Starting Raw Data Download and Organization ---"
    
    # Configuration for raw file download
    local base_url="https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed"
    local histone_marks=(
        "DNase" "H2A.Z" "H3K27ac" "H3K27me3" "H3K4me1" "H3K4me2"
        "H3K4me3" "H3K79me2" "H3K9ac" "H3K9me3" "H4K20me1"
    )
    local cell_type_count=129

    # --- Verification Step ---
    echo "Verifying remote file availability using E001 as a test case..."
    for mark in "${histone_marks[@]}"; do
        local filename="E001-${mark}.imputed.pval.signal.bigwig"
        local url="$base_url/$mark/$filename"
        if ! wget --spider -q "$url"; then
            echo "Error: Remote file not found at $url"
            echo "Aborting setup to prevent partial downloads."
            return 1 # Using return to exit the function with an error code
        fi
    done
    echo "Remote file verification successful."

    # --- Download Step ---
    mkdir -p "$RAW_DIR"
    echo "Creating subdirectories for each histone mark in $RAW_DIR..."
    for mark in "${histone_marks[@]}"; do
        mkdir -p "$RAW_DIR/$mark"
    done

    for i in $(seq 1 $cell_type_count); do
        local eid_num=$(printf "%03d" $i)
        local eid="E${eid_num}"
        echo "Processing data for cell type: $eid"
        
        for mark in "${histone_marks[@]}"; do
            local filename="${eid}-${mark}.imputed.pval.signal.bigwig"
            local target_file="$RAW_DIR/$mark/$filename"
            
            if [ -f "$target_file" ]; then
                echo "  - Found ${filename}, skipping."
            else
                local url="$base_url/$mark/$filename"
                echo "  - Downloading ${filename}..."
                if wget --spider -q "$url"; then
                    wget -q -O "$target_file" "$url"
                else
                    echo "    ...Warning: File not found on server at $url. Skipping."
                fi
            fi
        done
    done
    echo "--- Raw Data Download Finished ---"
    return 0
}

annotate_background() {
    local source_filepath="$1"
    local filename=$(basename "$source_filepath")
    local base_name="${filename%.tsv}"
    
    echo "Annotating background: $base_name with histone data..."

    time python3 "$ANNOTATE_SCRIPT_PY" \
        "$source_filepath" \
        "$BACKGROUND_DIR/${base_name}_annotated.parquet" \
        "$CACHED_DIR/${base_name}_cached.parquet" \
        "$RAW_DIR"
}

# --- Main Logic ---
echo "Starting Histone module setup..."

if [[ -z "$COMMON_BACKGROUND" || ! -d "$COMMON_BACKGROUND" ]]; then
    echo "Error: COMMON_BACKGROUND environment variable is not set or not a valid directory."
    exit 1
fi

shopt -s nullglob
source_files=("$COMMON_BACKGROUND"/*.tsv)
shopt -u nullglob

if [ ${#source_files[@]} -eq 0 ]; then
    echo "No .tsv files found in $COMMON_BACKGROUND. Nothing to set up."
    exit 1
fi

for filepath in "${source_files[@]}"; do
    filename=$(basename "$filepath")
    base_name="${filename%.tsv}"
    final_output_file="$BACKGROUND_DIR/${base_name}_annotated.parquet"

    echo "Processing manifest: $filename"

    if [ -f "$final_output_file" ]; then
        echo "  -> Annotated file already exists. Skipping."
        continue
    fi

    if [ -f "$PRECOMPUTED_URL_FILE" ]; then
        precomputed_base_url=$(cat "$PRECOMPUTED_URL_FILE")
        precomputed_url="${precomputed_base_url}${base_name}_annotated.parquet"
        echo "  -> Attempting to download pre-computed file from -> $precomputed_url <-"
        if wget -q --spider "$precomputed_url"; then
            wget -q -O "$final_output_file" "$precomputed_url"
            echo "  -> Pre-computed file downloaded successfully."
            continue
        else
            echo "  -> Pre-computed file not found. Falling back to local processing."
        fi
    fi

    echo "  -> Processing from raw data..."
    if [ ! -d "$RAW_DIR" ] || [ -z "$(ls -A "$RAW_DIR"/*/* 2>/dev/null)" ]; then
        echo "  -> Raw data not found locally. Attempting to download..."
        download_and_organize_raw_data || { echo "  -> CRITICAL: Failed to download raw data. Aborting setup for $filename."; continue; }
    else
        echo "  -> Raw data found locally."
    fi
    
    annotate_background "$filepath"
done

if [ -d "$RAW_DIR" ]; then
    echo "Cleaning up raw data directory..."
    rm -rf "$RAW_DIR"
fi

echo "Histone module setup complete."
