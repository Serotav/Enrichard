SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
PYTHON_GROUP_SCRIPT="$SCRIPT_DIR/analyze_groups.py"
USER_DIR=$1
METHOD=$2

python3 $PYTHON_GROUP_SCRIPT  --analysis-dir $USER_DIR 