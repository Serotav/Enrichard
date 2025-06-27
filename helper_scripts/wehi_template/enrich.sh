SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"

USER_DIR="$1"
OUPUT_DIR="$2"
BACKGROUND_NAME="$3"
P_VALUE="$4" 
CORRECTION="$5"

time python3 $SCRIPT_DIR/do_enrich.py \
    --user_dir "$USER_DIR" \
    --background_name "$BACKGROUND_NAME" \
    --p_value "$P_VALUE" \
    --correction "$CORRECTION" \
    --output_folder "$OUPUT_DIR" \