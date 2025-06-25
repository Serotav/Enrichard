SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"

USER_DIR="$1"
BACKGROUND_NAME="$2"
P_VALUE="$3" 
CORRECTION="$4"

time python3 $SCRIPT_DIR/do_enrich.py \
    --user_dir "$USER_DIR" \
    --background_name "$BACKGROUND_NAME" \
    --p_value "$P_VALUE" \
    --correction "$CORRECTION" 