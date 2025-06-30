SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"

TRAIT_TO_DOWNLOAD=$(cat "$SCRIPT_DIR/source_url.txt")

BACKGROUND_DIR="$SCRIPT_DIR/background"
RAW_BACKGROUND_DIR="$SCRIPT_DIR/raw_background"
ANNOTATE_SCRIPT_PY="$SCRIPT_DIR/annotate.py"
RDATA_DIR="$SCRIPT_DIR/rdata"
ANNOTATE_SCRIPT_R="$SCRIPT_DIR/GL2EntrezID.r"

mkdir -p "$RAW_BACKGROUND_DIR" "$BACKGROUND_DIR" "$RDATA_DIR" 

annotate_background() {
    local source_filepath="$1"
    local filename=$(basename "$source_filepath")
    
    cp "$source_filepath" "$RAW_BACKGROUND_DIR/$filename"

    Rscript $ANNOTATE_SCRIPT_R "$RAW_BACKGROUND_DIR/${filename%.gz}" "$RAW_BACKGROUND_DIR/${filename%.tsv.gz}_EnterezId.tsv" 2> /dev/null
    time python3 $ANNOTATE_SCRIPT_PY \
    "$RAW_BACKGROUND_DIR/${filename%.tsv.gz}_EnterezId.tsv" \
    "$BACKGROUND_DIR/${filename%.tsv.gz}_annotated.parquet"\
    $RDATA_DIR/* 

}

# Check if BACKGROUND_DIR is empty if not then skip this
if [ -d "$BACKGROUND_DIR" ] && [ -n "$(find "$BACKGROUND_DIR" -mindepth 1 -print -quit)" ]; then
    echo $SCRIPT_DIR "Background directory is already populated. Exiting."
    exit 0
fi

# Check that COMMON_BACKGROUND directory exists
if [[ -z "$COMMON_BACKGROUND" || ! -d "$COMMON_BACKGROUND" ]]; then
    echo "Error: COMMON_BACKGROUND environment variable is not set or not a valid directory."
    exit 1
fi

# Get all TSV files in the COMMON_BACKGROUND directory
shopt -s nullglob
source_files=("$COMMON_BACKGROUND"/*.tsv)
shopt -u nullglob 

# Download MSigDB Data Hallmark
wget -q $TRAIT_TO_DOWNLOAD -P "$RDATA_DIR"

# Parse Background Files
for filepath in "${source_files[@]}"; do
    if [ "$RUN_MODE" = "server" ]; then
        annotate_background "$filepath" &
    else
        annotate_background "$filepath"
    fi
done

wait

rm $RAW_BACKGROUND_DIR/* $RDATA_DIR/*
