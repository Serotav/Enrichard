#!/bin/bash
set -e
SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"

CHROMATINE_URL_FILE="$SCRIPT_DIR/source_url.txt"
CHROMATINE_RAW_DIR="$SCRIPT_DIR/raw_chromatine"
METADATA_FILE="$SCRIPT_DIR/EID_metadata.tab"

BACKGROUND_DIR="$SCRIPT_DIR/background"
CACHED_DIR="$SCRIPT_DIR/cached"
ANNOTATE_SCRIPT_PY="$SCRIPT_DIR/annotate.py"

# Create necessary directories
mkdir -p "$BACKGROUND_DIR" "$CACHED_DIR" "$CHROMATINE_RAW_DIR"

# Exit if background files are already generated
if [ -n "$(find "$BACKGROUND_DIR" -mindepth 1 -print -quit 2>/dev/null)" ]; then
    echo "Chromatine State background files already exist. Skipping setup."
    exit 0
fi

# Download and extract data if raw data directory is empty
if [ ! -n "$(find "$CHROMATINE_RAW_DIR" -mindepth 1 -print -quit 2>/dev/null)" ]; then
    if [ ! -f "$CHROMATINE_URL_FILE" ]; then
        echo "Error: source_url.txt not found!"
        exit 1
    fi
    CHROMATINE_URL=$(cat "$CHROMATINE_URL_FILE")
    echo "Downloading and extracting chromatin data from $CHROMATINE_URL..."
    wget -qO- "$CHROMATINE_URL" | tar -xzf - -C "$CHROMATINE_RAW_DIR"
else
    echo "Raw chromatin data already exists. Skipping download."
fi

echo "Downloading finishing chromatin state setup..."

# Check that COMMON_BACKGROUND environment variable is set
if [[ -z "$COMMON_BACKGROUND" || ! -d "$COMMON_BACKGROUND" ]]; then
    echo "Error: COMMON_BACKGROUND environment variable is not set or not a valid directory."
    exit 1
fi

# --- Annotation Function ---
annotate_background() {
    local source_filepath="$1"
    local filename=$(basename "$source_filepath")
    local base_name="${filename%.tsv}"
    
    echo "Annotating background: $base_name with chromatin states..."

    time python3 "$ANNOTATE_SCRIPT_PY" \
        "$source_filepath" \
        "$BACKGROUND_DIR/${base_name}_annotated.parquet" \
        "$CACHED_DIR/${base_name}_cached.parquet" \
        "$CHROMATINE_RAW_DIR" \
        "$METADATA_FILE"
}

# --- Main Loop ---
shopt -s nullglob
source_files=("$COMMON_BACKGROUND"/*.tsv)
shopt -u nullglob

if [ ${#source_files[@]} -eq 0 ]; then
    echo "No .tsv files found in $COMMON_BACKGROUND. Nothing to annotate."
    exit 1
fi

for filepath in "${source_files[@]}"; do
    # Use & to run in parallel if in server mode
    if [ "$RUN_MODE" = "server" ]; then
        annotate_background "$filepath" &
    else
        annotate_background "$filepath"
    fi
done

wait

# Clean up raw data to save space
rm -rf "$CHROMATINE_RAW_DIR"

echo "Chromatine State module setup complete."
