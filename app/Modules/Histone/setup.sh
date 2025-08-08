#!/bin/bash
set -e
SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
BACKGROUND_DIR="$SCRIPT_DIR/histone_annotated"

# Check if background directory exists
if [ ! -d "$BACKGROUND_DIR" ]; then
	echo "Error: Background directory '$BACKGROUND_DIR' does not exist."
	echo "READ the README.md of this module, make sure git lfs is setup on your machine"
	exit 1
fi


