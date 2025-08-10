import streamlit as st
import uuid
import shutil
import pathlib
import zipfile
import random
import polars as pl

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

def handle_zip_extraction(zip_file, extract_to_dir: pathlib.Path) -> None:
    """Extracts a zip file to a specified directory."""
    try:
        with zipfile.ZipFile(zip_file, 'r') as zip_ref:
            zip_ref.extractall(extract_to_dir)
    except Exception as e:
        st.error(f"Failed to extract zip file: {e}")
        st.stop()

def organize_extracted_samples(source_dir: pathlib.Path) -> None:
    """
    Recursively finds all sample files, moves them to a structured top-level directory,
    and renames them to the standard 'sample.txt' name.
    """
    
    # Find all files, excluding system files like .DS_Store
    all_files = [p for p in source_dir.rglob('*') if p.is_file() and not p.name.startswith('.')]
    
    if not all_files:
        st.warning("No sample files found in the extracted archive.")
        st.stop()

    for old_path in all_files:
        # Create a new directory named after the original file stem
        new_sample_dir = source_dir / old_path.stem
        new_sample_dir.mkdir(exist_ok=True)
        
        # Define the new path for the sample file
        new_path = new_sample_dir / USER_SAMPLE_NAME
        
        # Move and rename the file
        shutil.move(str(old_path), str(new_path))

    # Clean up any empty directories left behind after moving files
    for path in source_dir.iterdir():
        if path.is_dir() and not any(path.iterdir()):
            shutil.rmtree(path)
            


import polars as pl

def create_control_samples(real_samples_dir: pathlib.Path, background_path: pathlib.Path, control_samples_dir: pathlib.Path) -> None:
    """
    Generates a control sample for each real sample by random sampling from the background's 'Probe_ID' column.
    Assumes real_samples_dir contains one directory per sample.
    """
    try:
        control_samples_dir.mkdir(parents=True, exist_ok=True)
        
        # Use Polars to read the TSV and get the 'Probe_ID' column
        background_df = pl.read_csv(background_path, separator='\t', null_values='NA')
        if "Probe_ID" not in background_df.columns:
            st.error(f"The background file '{background_path.name}' must contain a 'Probe_ID' column.")
            st.stop()
        
        background_probes = background_df["Probe_ID"].to_list()

        if not background_probes:
            st.error("Background file does not contain any probes in the 'Probe_ID' column.")
            st.stop()

        # Iterate through the structured directories in real_samples_dir
        real_sample_dirs = [d for d in real_samples_dir.iterdir() if d.is_dir()]

        with st.spinner(f"Generating control samples for {len(real_sample_dirs)} real samples..."):
            for sample_dir in real_sample_dirs:
                real_sample_file = sample_dir / USER_SAMPLE_NAME
                if not real_sample_file.exists():
                    continue

                with open(real_sample_file, 'r') as f:
                    num_probes = sum(1 for line in f if line.strip())
                
                # Ensure we don't request more probes than available
                if num_probes > len(background_probes):
                    st.error(f"Sample '{sample_dir.name}' contains more probes ({num_probes}) than the background ({len(background_probes)}). Cannot create control.")
                    st.stop()

                control_probes = random.sample(background_probes, k=num_probes)
                
                # Create a corresponding directory structure for the control sample
                control_sample_dir = control_samples_dir / sample_dir.name
                control_sample_dir.mkdir(exist_ok=True)
                
                control_file_path = control_sample_dir / USER_SAMPLE_NAME
                with open(control_file_path, 'w') as f:
                    f.write('\n'.join(control_probes))
        
        st.success("Control samples generated successfully.")

    except FileNotFoundError:
        st.error(f"Background file not found at: {background_path}")
        st.stop()
    except pl.exceptions.NoDataError:
        st.error(f"The background file '{background_path.name}' appears to be empty.")
        st.stop()
    except ValueError as e:
        st.error(f"Error during random sampling: {e}.")
        st.stop()
    except Exception as e:
        st.error(f"An unexpected error occurred while creating control samples: {e}")
        st.stop()