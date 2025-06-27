import streamlit as st
import uuid
import shutil
import pathlib

from .config import *

def initialize_session()-> None:
    """
    Initializes the Streamlit session state.
    - Creates a unique user ID if one doesn't exist.
    - Sets up the user-specific directory path.
    """
    if 'user_id' not in st.session_state:
        st.session_state.user_id = str(uuid.uuid4())
    
    if 'user_dir' not in st.session_state:
        st.session_state.user_dir = USER_DATA_ROOT / st.session_state.user_id
        st.session_state.user_dir.mkdir(parents=True, exist_ok=True)

def get_background_options() -> list[str]:
    """
    Scans the COMMON_BACKGROUND_ROOT directory for available background files.
    Returns a list of unique names, split at the first dot.
    For example, 'HM450.hg38.manifest.tsv' becomes 'HM450'.
    """
    if not COMMON_BACKGROUND_ROOT.is_dir():
        st.error(f"Configuration Error: Background directory '{COMMON_BACKGROUND_ROOT}' not found.")
        return []
    
    options = [] 
    
    for file_path in COMMON_BACKGROUND_ROOT.glob('*.tsv'):
        first_part = file_path.name.split('.', 1)[0]
        options.append(first_part)

    return sorted(options)

def get_example_files() -> dict:
    """
    Scans the TEST_SAMPLES_DIR for available example files.
    Returns a dictionary mapping a display name to its file path.
    """
    if not TEST_SAMPLES_DIR.is_dir():
        return {} 
    
    example_files = {
        file_path.stem: file_path 
        for file_path in TEST_SAMPLES_DIR.glob('*')
    }
    return example_files

def save_uploaded_file(uploaded_file, destination_path: pathlib.Path) -> None:
    """Saves an uploaded file to a specified destination."""
    try:
        with open(destination_path, "wb") as f:
            f.write(uploaded_file.getvalue())
    except Exception as e:
        st.error(f"Failed to save file '{uploaded_file.name}': {e}")
        st.stop()

def copy_example_file(source_path: pathlib.Path, destination_path: pathlib.Path) -> None:
    """Copies an example file to the user's session directory."""
    try:
        shutil.copy(source_path, destination_path)
    except Exception as e:
        st.error(f"Failed to copy example file '{source_path.name}': {e}")
        st.stop()