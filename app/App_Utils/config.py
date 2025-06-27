import os
import pathlib

APP_DIR = pathlib.Path(__file__).parent.parent

print("App directory:", APP_DIR)
TEST_SAMPLES_DIR = APP_DIR / "Test_samples"
USER_DATA_ROOT = APP_DIR / "User_Data"
COMMON_BACKGROUND_ROOT = pathlib.Path(os.getenv("COMMON_BACKGROUND"))
MASTER_SCRIPT_PATH = APP_DIR / "master_enrich.sh"

USER_CUSTOM_BACKGROUND_NAME = os.getenv("USER_CUSTOM_BACKGROUND_NAME")
USER_SAMPLE_NAME = os.getenv("USER_SAMPLE_NAME")
USER_MODULE_OUTPUT = os.getenv("USER_MODULE_OUTPUT")

TOW_SAMPLE_COMPARISON_NAME_1 = "Sample_A"
TOW_SAMPLE_COMPARISON_NAME_2 = "Sample_B"