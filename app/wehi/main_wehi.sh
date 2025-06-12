#!/bin/bash

BASE_DIR=$(dirname "$0")

# Define relative paths
RDATA_DIR="$BASE_DIR/rdata"
BACKGROUND_DIR="$BASE_DIR/background"
SAMPLE_DIR="$BASE_DIR/sample"
PARSED_DIR="$BASE_DIR/Parsed"

# Create directories if they don't exist
mkdir -p "$RDATA_DIR" "$BACKGROUND_DIR" "$SAMPLE_DIR" "$PARSED_DIR"

ANNOTATE_SCRIPT_R="$BASE_DIR/look_up_table_create.r"
ANNOTATE_SCRIPT_PY="$BASE_DIR/annotate.py"
ENRICH_SCRIPT="$BASE_DIR/do_enrich.py"

# --- Input/Output Files --
SAMPLE_FILE="$SAMPLE_DIR/uploaded_sample.csv"
BACKGROUND_FILE="$BACKGROUND_DIR/MSA_parsed.csv"
OUTPUT_FILE="$PARSED_DIR/significant_results.csv"
MAPPING_FILE="$BACKGROUND_DIR/symbols_to_entrezid.txt"

# --- Download MSigDB Data ---
if [ ! -f "$RDATA_DIR/human_H_v5p2.rdata" ]; then
    echo "Downloading MSigDB .rdata files..."
    wget https://bioinf.wehi.edu.au/MSigDB/v7.1/Hs.h.all.v7.1.entrez.rds -P "$RDATA_DIR"
    wget https://bioinf.wehi.edu.au/MSigDB/v7.1/Hs.c{1..7}.all.v7.1.entrez.rds -P "$RDATA_DIR"
fi

# --- Create Gene ID Mapping ---
if [ ! -f "$MAPPING_FILE" ]; then
    Rscript "$ANNOTATE_SCRIPT_R"
    if [ $? -ne 0 ]; then
        echo "Error running R script ($ANNOTATE_SCRIPT_R). Exiting."
        exit 1
    fi
fi

# --- Create Background Annotation File ---
if [ ! -f "$BACKGROUND_FILE" ]; then
    python3 "$ANNOTATE_SCRIPT_PY"
    if [ $? -ne 0 ]; then
        echo "Error running Python script ($ANNOTATE_SCRIPT_PY). Exiting."
        exit 1
    fi
fi

# --- Run Enrichment Analysis ---
echo "Running enrichment analysis..."
if [ ! -f "$SAMPLE_FILE" ]; then
    echo "Error: Sample file ($SAMPLE_FILE) not found. Did the upload succeed? Exiting."
    exit 1
fi

python3 "$ENRICH_SCRIPT" "$BACKGROUND_FILE" "$SAMPLE_FILE" "$OUTPUT_FILE" --p_value_threshold 0.05

if [ $? -ne 0 ]; then
    echo "Error running enrichment script ($ENRICH_SCRIPT). Check logs."
    exit 1
fi

exit 0