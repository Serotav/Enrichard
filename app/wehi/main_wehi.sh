#!/bin/bash

BASE_DIR=$(dirname "$0")
BACKGROUND_CHOICE=$1

# Define relative paths
RDATA_DIR="$BASE_DIR/rdata"
BACKGROUND_DIR="$BASE_DIR/background"
SAMPLE_DIR="$BASE_DIR/sample"
PARSED_DIR="$BASE_DIR/Parsed"

# Create directories if they don't exist
mkdir -p  "$SAMPLE_DIR" "$PARSED_DIR"

ANNOTATE_SCRIPT_PY="$BASE_DIR/annotate.py"
ENRICH_SCRIPT="$BASE_DIR/do_enrich.py"

# --- Input/Output Files --
SAMPLE_FILE="$SAMPLE_DIR/uploaded_sample.csv"
OUTPUT_FILE="$PARSED_DIR/significant_results.csv"

# --- Run Enrichment Analysis ---
echo "Running enrichment analysis..."
if [ ! -f "$SAMPLE_FILE" ]; then
    echo "Error: Sample file ($SAMPLE_FILE) not found. Did the upload succeed? Exiting."
    exit 1
fi

python3 "$ENRICH_SCRIPT" "${BACKGROUND_CHOICE}" "$SAMPLE_FILE" "$OUTPUT_FILE" --p_value_threshold 0.05

if [ $? -ne 0 ]; then
    echo "Error running enrichment script ($ENRICH_SCRIPT). Check logs."
    exit 1
fi

exit 0