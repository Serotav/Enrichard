import streamlit as st
import pandas as pd
import subprocess
import os

# --- Configuration ---
# Define paths relative to this app.py file
APP_DIR = os.path.dirname(__file__)
WEHI_DIR = os.path.join(APP_DIR, "wehi")
SAMPLE_DIR = os.path.join(WEHI_DIR, "sample")
PARSED_DIR = os.path.join(WEHI_DIR, "Parsed")
MAIN_WEHI_SCRIPT = "main_wehi.sh"

UPLOADED_SAMPLE_PATH = os.path.join(SAMPLE_DIR, "uploaded_sample.csv")
UPLOADED_BACKGROUND_PATH = os.path.join(SAMPLE_DIR, "uploaded_background.csv")
RESULTS_PATH = os.path.join(PARSED_DIR, "significant_results.csv")
BACKGROUND_DIR = os.path.join(WEHI_DIR, "background")

os.makedirs(SAMPLE_DIR, exist_ok=True)
os.makedirs(PARSED_DIR, exist_ok=True)

# Scans a directory for background files and returns a dictionary mapping display names to full filenames.
def get_background_options(directory):
    options = {}
    if os.path.isdir(directory):
        for filename in os.listdir(directory):
            if filename.endswith(('.tsv', '.csv')): 
                display_name = filename.split('.')[0]
                options[display_name] = filename
    else:
        raise FileNotFoundError(f"Background directory '{directory}' does not exist or is not a directory.")
    return options


# --- Streamlit App ---
st.set_page_config(layout="wide")
st.title("CpG Site Enrichment Analysis")

st.markdown("""
Upload a file containing CpG probe IDs (one ID per line, no header)
to perform enrichment analysis using WEHI MSigDB gene sets.
""")

uploaded_sample_file = st.file_uploader("Choose a file with your CpG sites (.csv or .txt)", type=['csv', 'txt'], help="One CpG ID per line, no header.")

if uploaded_sample_file is not None:
    # --- Step 2: Select the background universe ---
    st.header("2. Select Background")
    
    wehi_backgrounds = get_background_options(BACKGROUND_DIR)

    selected_background = st.selectbox(
        "Choose the background set for enrichment testing:",
        options=list(wehi_backgrounds.keys()) + ['custom'],
        index=0
    )

    custom_background_file = None
    # Conditionally display the file uploader for a custom background
    if selected_background == 'custom':
        st.info("Please upload a .txt file containing the list of CpG sites for your custom background.")
        custom_background_file = st.file_uploader(
            "Upload your custom background file (.txt)", 
            type=['txt','csv'], 
            key='custom_bg_uploader' # Unique key for this uploader
        )
    
    st.markdown("---")
    
    # --- Step 3: Run the analysis ---
    st.header("3. Run Analysis")
    if st.button("Run Enrichment Analysis"):
        # --- Pre-run checks and file saving ---
        try:
            # Save the main sample file
            with open(UPLOADED_SAMPLE_PATH, "wb") as f:
                f.write(uploaded_sample_file.getvalue())
            st.success(f"Sample file '{uploaded_sample_file.name}' saved.")

            # Build the command to run the shell script
            command_to_run = ['bash', MAIN_WEHI_SCRIPT]

            # Handle the custom background case
            if selected_background == 'custom':
                if custom_background_file is None:
                    st.error("You selected 'custom' background but did not upload a file. Please upload a background file.")
                    st.stop() # Stop execution if the file is missing
                
                # Save the custom background file
                with open(UPLOADED_BACKGROUND_PATH, "wb") as f:
                    f.write(custom_background_file.getvalue())
                st.success(f"Custom background file '{custom_background_file.name}' saved.")
                
                # Add the path of the custom background file as an argument
                command_to_run.append(UPLOADED_BACKGROUND_PATH)

            else:
                # Use the selected background file from the predefined options
                command_to_run.append(os.path.join(BACKGROUND_DIR, wehi_backgrounds[selected_background]))
                st.success(f"Using background file: {selected_background}")

            # --- Execute the pipeline ---
            st.info(f"Starting analysis pipeline with background: '{selected_background}'")
            st.markdown("---")

            # Clean up previous results if they exist
            if os.path.exists(RESULTS_PATH):
                os.remove(RESULTS_PATH)

            # Run the main_wehi.sh script with the correct arguments
            with st.spinner(f"Running `main_wehi.sh`... This may take a few minutes."):
                process = subprocess.run(
                    command_to_run,  # Use the dynamically built command
                    capture_output=True,
                    text=True,
                    cwd=WEHI_DIR 
                )

            st.markdown("### Pipeline Output Log")
            st.code(f"Command executed: {' '.join(command_to_run)}\n\nReturn Code: {process.returncode}\n\nSTDOUT:\n{process.stdout}\n\nSTDERR:\n{process.stderr}", language='log')

            if process.returncode == 0:
                st.success("Pipeline completed successfully!")
                # Try to read and display results
                if os.path.exists(RESULTS_PATH):
                    try:
                        results_df = pd.read_csv(RESULTS_PATH)
                        st.markdown("### Significant Enrichment Results (P < 0.05)")
                        if not results_df.empty:
                             st.dataframe(results_df)
                        else:
                             st.info("No statistically significant enrichment found for the provided CpG sites.")
                    except pd.errors.EmptyDataError:
                         st.info("No statistically significant enrichment found (results file is empty).")
                    except Exception as e:
                        st.error(f"Error reading results file ({RESULTS_PATH}): {e}")
                else:
                     st.warning("Pipeline finished, but the results file was not found.")

            else:
                st.error("Pipeline execution failed. Check the log above for details.")

        except Exception as e:
            st.error(f"An error occurred during the process: {e}")
else:
    st.info("Please upload a file to begin.")

st.markdown("---")