#!/bin/bash
SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"

USER_DIR="$1"
OUPUT_DIR="$2"
BACKGROUND_NAME="$3"
P_VALUE="$4" 
CORRECTION="$5"

BACKGROUND_DIR="$SCRIPT_DIR/histone_annotated"

# Check if background directory exists
if [ ! -d "$BACKGROUND_DIR" ]; then
	echo "Error: Background directory '$BACKGROUND_DIR' does not exist."
	echo "READ the README.md of this module, make sure git lfs is setup on your machine"
	exit 1
fi


time python3 "$SCRIPT_DIR/do_enrich.py" \
    --user_dir "$USER_DIR" \
    --background_name "$BACKGROUND_NAME" \
    --p_value "$P_VALUE" \
    --correction "$CORRECTION" \
    --output_folder "$OUPUT_DIR"
