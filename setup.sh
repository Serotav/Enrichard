#!/bin/bash

BASE_DIR=$(dirname "$0")/app/wehi
BACKGROUND_DIR="$BASE_DIR/background"
RAW_BACKGROUND_DIR="$BASE_DIR/raw_background"
ANNOTATE_SCRIPT_PY="$BASE_DIR/annotate.py"
RDATA_DIR="$BASE_DIR/rdata"
ANNOTATE_SCRIPT_R="$BASE_DIR/GL2EntrezID.r"

mkdir -p "$RAW_BACKGROUND_DIR" "$BACKGROUND_DIR" "$RDATA_DIR" 

parse() {
    local url="$1"
    local filename=$(basename "$url")
    
    wget -q "$url" -O "$RAW_BACKGROUND_DIR/$filename"
    gzip -d "$RAW_BACKGROUND_DIR/$filename"
    # BRO R NEEDS TO SHUT THE FUCK UP AND NOT MESS MY LOGS
    Rscript $ANNOTATE_SCRIPT_R "$RAW_BACKGROUND_DIR/${filename%.gz}" "$RAW_BACKGROUND_DIR/${filename%.tsv.gz}_EnterezId.tsv" 2> /dev/null
    python3 $ANNOTATE_SCRIPT_PY "$RAW_BACKGROUND_DIR/${filename%.tsv.gz}_EnterezId.tsv" \
    "$BACKGROUND_DIR/${filename%.tsv.gz}_annotated.tsv"\
    $RDATA_DIR
}

URLS=(
    "https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/MSA/MSA.hg38.manifest.tsv.gz"
    "https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/EPICv2/EPICv2.hg38.manifest.tsv.gz"
    "https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/EPIC+/EPIC+.hg38.manifest.tsv.gz"
    "https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/EPIC/EPIC.hg38.manifest.tsv.gz"
    "https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/HM450/HM450.hg38.manifest.tsv.gz"
)

# --- Download Background ---
if [ -z "$(ls -A "$BACKGROUND_DIR" 2>/dev/null)" ]; then
    for url in "${URLS[@]}"; do
        parse "$url" &
    done
fi

# --- Download MSigDB Data ---
if [ ! -f "$RDATA_DIR/Hs.h.all.v7.1.entrez.rds" ]; then
    echo "Downloading MSigDB .rdata files..."
    wget https://bioinf.wehi.edu.au/MSigDB/v7.1/Hs.h.all.v7.1.entrez.rds -P "$RDATA_DIR"
    wget https://bioinf.wehi.edu.au/MSigDB/v7.1/Hs.c{1..7}.all.v7.1.entrez.rds -P "$RDATA_DIR"
fi



wait
rm $RAW_BACKGROUND_DIR/*
