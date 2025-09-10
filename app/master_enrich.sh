SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
MODULES_DIR="$SCRIPT_DIR/Modules"

USER_DIR="$1"
OUPUT_DIR="$2"
BACKGROUND_NAME="$3"
P_VALUE="$4" 
CORRECTION="$5"

echo "Running master enrichment script with parameters:" $USER_DIR >&2

# Call enrioch script for each module
for module_dir in $MODULES_DIR/*/; do
    module_name=$(basename "$module_dir")
    echo "--- Starting module: $module_name ---" >&2

    if [ -f "$module_dir/enrich.sh" ]; then
        if [ "$RUN_MODE" = "server" ]; then
            (
                time bash "$module_dir/enrich.sh" $USER_DIR $OUPUT_DIR $BACKGROUND_NAME $P_VALUE $CORRECTION
                echo "--- Finished module: $module_name ---"
                echo 
            ) >&2 &
        else
            time bash "$module_dir/enrich.sh" $USER_DIR $OUPUT_DIR $BACKGROUND_NAME $P_VALUE $CORRECTION
        fi
    fi
done
wait 