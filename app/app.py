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
RESULTS_PATH = os.path.join(PARSED_DIR, "significant_results.csv")

os.makedirs(SAMPLE_DIR, exist_ok=True)
os.makedirs(PARSED_DIR, exist_ok=True)

# --- Streamlit App ---
st.set_page_config(layout="wide")
st.title("CpG Site Enrichment Analysis")

st.markdown("""
Upload a file containing CpG probe IDs (one ID per line, no header)
to perform enrichment analysis using WEHI MSigDB gene sets.
""")

uploaded_file = st.file_uploader("Choose a file (.csv or .txt)", type=['csv', 'txt'])

if uploaded_file is not None:
    try:
        with open(UPLOADED_SAMPLE_PATH, "wb") as f:
            f.write(uploaded_file.getvalue())
        st.success(f"File '{uploaded_file.name}' saved successfully.")

        if st.button("Run Enrichment Analysis"):
            st.info("Starting analysis pipeline...")
            st.markdown("---") 

            if os.path.exists(RESULTS_PATH):
                os.remove(RESULTS_PATH)

            # Run the main_wehi.sh script
            with st.spinner("Running `main_wehi.sh` ..."):
                process = subprocess.run(
                    ['bash', MAIN_WEHI_SCRIPT],
                    capture_output=True,
                    text=True,
                    cwd=WEHI_DIR 
                )

            st.markdown("### Pipeline Output Log")
            st.code(f"Return Code: {process.returncode}\n\nSTDOUT:\n{process.stdout}\n\nSTDERR:\n{process.stderr}", language='log')

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
                     st.warning("Pipeline finished, but the results file was not found.") # Should not happen if script logic is correct

            else:
                st.error("Pipeline execution failed. Check the log above for details.")

    except Exception as e:
        st.error(f"An error occurred processing the uploaded file: {e}")
else:
    st.info("Please upload a file to begin.")


st.markdown("---")