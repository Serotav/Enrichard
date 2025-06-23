SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
MODULES_DIR="$SCRIPT_DIR/Modules"

USER_DIR="$1"
BACKGROUND_NAME="$2"
P_VALUE="$3" 

# Call enrioch script for each module
for module_dir in $MODULES_DIR/*/; do
    echo "Running enrich for $(basename "$module_dir")"
    if [ -f "$module_dir/setup.sh" ]; then
        bash "$module_dir/enrich.sh" $USER_DIR $BACKGROUND_NAME $P_VALUE # &
    fi
done

wait 