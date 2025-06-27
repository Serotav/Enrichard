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

def top_style():
    st.markdown(
        """
        <style>
            /* The banner styling */
            .fixed-header {
                position: fixed;
                top: 0;
                left: 0;
                width: 100%;
                background-color: #ecebe3;  
                border-bottom: 1px solid #d3d2ca; 
                color: #3d3a2a;             
                padding: 10px 0;
                z-index: 999;
                
                /* Flexbox for alignment */
                display: flex;
                justify-content: center;
                align-items: center;
            }
            
            /* Styling for the link within the banner */
            .fixed-header a {
                color: #3d3a2a;  /
                text-decoration: none;
                font-weight: 500;
            }

            .fixed-header a:hover {
                color: #bb5a38; 
                text-decoration: none;
            }
            
            /* Styling for the inline SVG GitHub logo */
            .github-logo {
                height: 22px;       
                width: 22px;
                margin-right: 10px; 
                fill: currentColor; 
            }

            /* Hide the default Streamlit menu */
            [data-testid="stToolbar"] {
                display: none !important;
            }

            /* Push the main app content down to avoid overlap */
            .main .block-container {
                padding-top: 5rem;
            }
        </style>
        """,
        unsafe_allow_html=True
    )
    
    st.markdown(
        """
        <div class="fixed-header">
            <a href="https://github.com/Serotav/Enrichard" target="_blank">Check out the source on GitHub</a>
        </div>
        """,
        unsafe_allow_html=True
    )