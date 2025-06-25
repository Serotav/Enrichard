SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
MODULES_DIR="$SCRIPT_DIR/Modules"

USER_DIR="$1"
BACKGROUND_NAME="$2"
P_VALUE="$3" 
CORRECTION="$4"

# Call enrioch script for each module
for module_dir in $MODULES_DIR/*/; do
    echo "Running enrich for $(basename "$module_dir")"
    if [ -f "$module_dir/setup.sh" ]; then
        if [ "$RUN_MODE" = "server" ]; then
            bash "$module_dir/enrich.sh" $USER_DIR $BACKGROUND_NAME $P_VALUE $CORRECTION &
        else
            bash "$module_dir/enrich.sh" $USER_DIR $BACKGROUND_NAME $P_VALUE $CORRECTION
        fi
    fi
    break
done

wait 