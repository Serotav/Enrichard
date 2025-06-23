#!/bin/bash

SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
MODULES_DIR="$SCRIPT_DIR/Modules"
USER_DATA="$SCRIPT_DIR/User_Data"

rm -rf "$USER_DATA"/*

URLS=(
    "https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/MSA/MSA.hg38.manifest.tsv.gz"
    #"https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/EPICv2/EPICv2.hg38.manifest.tsv.gz"
    #"https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/EPIC+/EPIC+.hg38.manifest.tsv.gz"
    #"https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/EPIC/EPIC.hg38.manifest.tsv.gz"
    "https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/HM450/HM450.hg38.manifest.tsv.gz"
)

# This is set in docker env
mkdir -p "$COMMON_BACKGROUND"

# Check if COMMON_BACKGROUND is empty if not then skip this
if [ -d "$COMMON_BACKGROUND" ] && [ -n "$(find "$COMMON_BACKGROUND" -mindepth 1 -print -quit)" ]; then
    echo $SCRIPT_DIR "Background directory is already populated. Exiting."
    exit 0
fi

# Download and decompress background files
for url in "${URLS[@]}"; do
    filename=$(basename "$url")
    wget -q "$url" -O "$COMMON_BACKGROUND/$filename" 
    gzip -d "$COMMON_BACKGROUND/$filename"
done

# Call setup scripts in each module directory
for module_dir in $MODULES_DIR/*/; do
    echo "Running setup for $(basename "$module_dir")"
    if [ -f "$module_dir/setup.sh" ]; then
        bash "$module_dir/setup.sh" # &
    fi
done

wait