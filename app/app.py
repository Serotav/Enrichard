import streamlit as st
import pandas as pd
import subprocess
import os
import pathlib
import uuid 
import shutil 
import altair as alt
import polars as pl
import numpy as np
from statsmodels.stats.multitest import multitest_methods_names

APP_DIR = pathlib.Path(__file__).parent
TEST_SAMPLES_DIR = APP_DIR / "Test_samples"
USER_DATA_ROOT = APP_DIR / "User_Data"
COMMON_BACKGROUND_ROOT = pathlib.Path(os.getenv("COMMON_BACKGROUND"))
MASTER_SCRIPT_PATH = APP_DIR / "master_enrich.sh"

USER_CUSTOM_BACKGROUND_NAME = os.getenv("USER_CUSTOM_BACKGROUND_NAME")
USER_SAMPLE_NAME = os.getenv("USER_SAMPLE_NAME")
USER_MODULE_OUTPUT = os.getenv("USER_MODULE_OUTPUT")

def initialize_session():
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

def create_dot_plot(df: pd.DataFrame):
    """
    Creates an Altair dot plot for enrichment analysis results.
    """
    # Only show the top N results to keep the chart readable
    # We can do .head() becouse they are already sorted by P-Value in the DataFrame
    df_to_plot = df.head(20).copy()

    min_p_val = df_to_plot['P-Value'].min()
    max_p_val = df_to_plot['P-Value'].max()

    # Create a list of 5 evenly-spaced values for the legend, 
    # ensuring the min and max from the data are included.
    legend_p_values = np.linspace(min_p_val, max_p_val, 5).tolist()

    # We need to sort the traits for the y-axis based on P-Value for a clean look
    sort_order = df_to_plot.sort_values("P-Value")["Trait"].tolist()

    chart = alt.Chart(df_to_plot).mark_circle().encode(
        y=alt.Y('Trait:N', sort=sort_order, title="Enriched Trait"),
        x=alt.X('Fold-Change:Q', title="Fold Change", scale=alt.Scale(zero=False)),
        color=alt.Color('P-Value:Q', 
                        scale=alt.Scale(scheme='yelloworangered', reverse=True), 
                        title="P-Value" ,legend=alt.Legend(format=".2e",values=legend_p_values)),
        
        size=alt.Size('a:Q', title="Count in Sample"),

        tooltip=[
            alt.Tooltip('Trait:N'),
            alt.Tooltip('P-Value:Q', format=".2e"), 
            alt.Tooltip('Fold-Change:Q', format=".2f"),
            alt.Tooltip('a:Q', title="Sample Hits")
        ]
    ).properties(
        title="Top Enriched Pathways/Traits"
    ).interactive()

    return chart

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

def run_enrichment_pipeline(user_dir: pathlib.Path, background: str, p_value: float, correction_method: str) -> None:
    """
    Executes the master enrichment shell script and displays its output.
    """
    command = ['bash', str(MASTER_SCRIPT_PATH), str(user_dir), str(background), str(p_value), str(correction_method)]
    
    st.markdown("---")

    with st.spinner("Running enrichment analysis. We'll be right back with results!"):
        process = subprocess.run(
            command,
            capture_output=True,
            text=True,
        )

    with st.expander("Show Pipeline Execution Log"):
        st.code(
            f"Command: {' '.join(command)}\n"
            f"Return Code: {process.returncode}\n\n"
            f"--- STDOUT ---\n{process.stdout}\n\n"
            f"--- STDERR ---\n{process.stderr}",
            language='log'
        )

    if process.returncode != 0:
        st.error("Pipeline execution failed. Please check the log above for details.")
        st.stop()
    else:
        st.success("Pipeline completed successfully!")

def render_methodology_explanation():
    """Renders the explanation of the Fisher's Exact Test methodology."""
    st.markdown("""
    For each trait, a **Fisher's Exact Test** was performed to determine if the trait is significantly enriched in your sample compared to the background universe. This test uses a 2x2 contingency table:
    """)
    
    col1, col2 = st.columns(2)
    with col1:
        st.markdown("""
        | | Has Trait | No Trait |
        |---|---|---|
        | **In Sample** | **a** | **b** |
        | **Not in Sample** | **c** | **d** |
        """)
    with col2:
        st.markdown("""
        - **a**: Probes in your sample that have the trait.
        - **b**: Probes in your sample that **do not** have the trait.
        - **c**: Probes in the background (but not your sample) that have the trait.
        - **d**: Probes in the background (but not your sample) that **do not** have the trait.
        """)
    
    st.markdown("The **P-Value** represents the probability of observing such an enrichment (or a greater one) by random chance. A lower P-Value indicates a more statistically significant enrichment.")
    st.markdown("The **Odds-Ratio** measures the strength of the association. It is calculated as:")
    st.latex(r'\text{Odds-Ratio} = \frac{a \times d}{b \times c}')
    st.markdown("An Odds-Ratio greater than 1 suggests that the trait is more likely to be found in your sample group than in the background.")

def display_single_module_results(module_name: str, file_path: pathlib.Path) -> bool:
    """
    Handles the logic for displaying the results (chart and data) for one module.
    Returns True if results were found and displayed, False otherwise.
    """
    st.markdown(f"#### Module: `{module_name}`")
    try:
        df = pl.read_csv(file_path, separator='\t').to_pandas()
        
        if df.empty:
            st.info("This module produced no significant results.")
            return False

        st.success(f"Found {len(df)} significant traits.")
        
        st.subheader("Visualization")
        dot_plot = create_dot_plot(df)
        st.altair_chart(dot_plot, use_container_width=True)
        
        with st.expander("Show Full Data Table"):
            float_cols = df.select_dtypes(include='float').columns
            format_dict = {col: '{:.2e}' for col in float_cols}
            st.dataframe(df.style.format(format_dict))
        
        return True

    except pl.exceptions.NoDataError:
        st.info("This module produced no significant results (empty file).")
        return False
    except Exception as e:
        st.error(f"Error reading or processing result file '{file_path.name}': {e}")
        return False

def display_results():
    """
    Scans for result files, creates tabs, and displays results and methodology.
    """
    module_output = st.session_state.user_dir / USER_MODULE_OUTPUT
    result_files = sorted(list(module_output.glob('*/*.tsv')))
    
    if not result_files:
        st.warning("Analysis complete, but no result files were found.")
        return

    st.header("Enrichment Results")

    with st.expander("How to Interpret These Results (Methodology)", expanded=False):
        render_methodology_explanation()
    
    st.markdown("---") 

    module_names = [path.parent.name for path in result_files]
    tabs = st.tabs(module_names)

    results_found_in_any_module = False
    for i, file_path in enumerate(result_files):
        with tabs[i]:
            if display_single_module_results(module_names[i], file_path):
                results_found_in_any_module = True

    if not results_found_in_any_module:
        st.info("The pipeline ran, but no modules found significant enrichment.")


def main():
    st.set_page_config(layout="wide")

    st.markdown("<h1 style='text-align: center;'>🚀 E N R I C H A R D 🚀</h1>", unsafe_allow_html=True)
    st.subheader("CpG Site Enrichment Analysis")
    initialize_session()

    st.markdown("Upload a file containing CpG probe IDs (one ID per line) to perform enrichment analysis.")
    
    input_method = st.radio(
        "Choose your input method:",
        ("Upload a File", "Use an Example"),
        horizontal=True,
        label_visibility="collapsed"
    )

    sample_source = None
    uploaded_or_example = None

    if input_method == "Upload a File":
        uploaded_sample_file = st.file_uploader(
            "Upload your list of CpG sites (one ID per line, no header)", type=['csv', 'txt']
        )
        if uploaded_sample_file:
            sample_source = uploaded_sample_file
            uploaded_or_example = 'uploaded'

    elif input_method == "Use an Example":
        example_files = get_example_files()
        if not example_files:
            st.info("No example files found.")
        else:
            options = ["- Select an example -"] + list(example_files.keys())
            selected_example_name = st.selectbox(
                "Choose an example file to run:",
                options=options
            )
            if selected_example_name != "- Select an example -":
                sample_source = example_files[selected_example_name]
                uploaded_or_example = 'example'
    

    if sample_source:
        # Configure Analysis Parameters 
        st.header("Configure Analysis Parameters")

        # Create two columns for the configuration options
        col1, col2, col3 = st.columns(3)

        custom_background_file = None  # Initialize to None

        #  Left Column: Background Selection 
        with col1:
            st.subheader("Background Universe")
            available_choices = get_background_options() + ['custom']
            selected_background_name = st.selectbox(
                label= "Choose the background set",
                options=available_choices,
            )
            
            if selected_background_name == 'custom':
                st.info("Your custom background should be a file with one CpG ID per line.")
                custom_background_file = st.file_uploader(
                    "Upload your custom background file", type=['csv', 'txt'], key='custom_bg'
                )
        
        # Center Column: p val
        with col2:
            st.subheader("Statistical Options")
            p_value_threshold = st.number_input(
                label="P-value threshold",
                min_value=0.0, max_value=1.0, value=0.05, step=0.01,
                format="%.2f", help="The significance threshold for the Fisher's Exact Test."
            )
        
        # Right Column: Correction Method
        with col3:
            st.subheader("Correction Method")
            # Get the list of methods and add 'none'
            correction_options = sorted(list(multitest_methods_names.values()) + ['none'])
            
            # Set a sensible default index ('fdr_bh' is highly recommended)
            default_index = correction_options.index('Bonferroni') 

            selected_correction_method = st.selectbox(
                label="Multiple Testing Correction",
                options=correction_options,
                index=default_index,
                help="Method to adjust p-values for multiple comparisons."
            )
        
        # Determine if the 'Run' button should be disabled
        is_run_disabled = (selected_background_name == 'custom' and custom_background_file is None)

        st.markdown("---") 

        if st.button("Run Enrichment Analysis", type="primary", disabled=is_run_disabled):
            # Clean up the directory from any previous run.
            user_dir = st.session_state.user_dir
            if user_dir.exists():
                shutil.rmtree(user_dir)
            user_dir.mkdir(exist_ok=True)
            
            # Handel both uploaded and example files
            user_sample_path = user_dir / USER_SAMPLE_NAME
            if uploaded_or_example == 'uploaded': 
                st.info(f"Using uploaded file: {sample_source.name}")
                save_uploaded_file(sample_source, user_sample_path)
            elif uploaded_or_example == 'example':
                st.info(f"Using example file: {sample_source.name}")
                copy_example_file(sample_source, user_sample_path)
            else:
                st.error("No valid sample file provided. Please upload or select an example file.")
                st.stop()

            background = selected_background_name 
            if selected_background_name == 'custom':
                background = USER_CUSTOM_BACKGROUND_NAME 
                user_background_path = user_dir / USER_CUSTOM_BACKGROUND_NAME
                save_uploaded_file(custom_background_file, user_background_path)

            # Run the enrichment analysis pipeline
            run_enrichment_pipeline(user_dir, background, p_value_threshold, selected_correction_method)
            display_results()
            
    else:
        st.info("Please upload a file to begin.")


if __name__ == "__main__":
    main()