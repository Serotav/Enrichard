import streamlit as st
import pandas as pd
import subprocess
import os
import pathlib
import uuid 
import shutil 
import altair as alt
import polars as pl

APP_DIR = pathlib.Path(__file__).parent
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

def create_dot_plot(df: pd.DataFrame):
    """
    Creates an Altair dot plot for enrichment analysis results.
    """
    # Only show the top N results to keep the chart readable
    # We can do .head() becouse they are already sorted by P-Value in the DataFrame
    df_to_plot = df.head(20).copy()

    # We need to sort the traits for the y-axis based on P-Value for a clean look
    sort_order = df_to_plot.sort_values("P-Value")["Trait"].tolist()

    chart = alt.Chart(df_to_plot).mark_circle().encode(
        y=alt.Y('Trait:N', sort=sort_order, title="Enriched Trait"),
        x=alt.X('Odds-Ratio:Q', title="Odds Ratio"),
        color=alt.Color('P-Value:Q', 
                        scale=alt.Scale(scheme='lightgreyred', reverse=True), 
                        title="P-Value"),
        
        size=alt.Size('a:Q', title="Count in Sample"),

        tooltip=[
            alt.Tooltip('Trait:N'),
            alt.Tooltip('P-Value:Q', format=".2e"), 
            alt.Tooltip('Odds-Ratio:Q', format=".2f"),
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

def run_enrichment_pipeline(user_dir: pathlib.Path, background: str, p_value: float) -> None:
    """
    Executes the master enrichment shell script and displays its output.
    """
    command = ['bash', str(MASTER_SCRIPT_PATH), str(user_dir), str(background), str(p_value)]
    
    st.info("Starting analysis pipeline...")
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
            st.dataframe(df)
        
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

    st.markdown("Upload a file containing CpG probe IDs (one ID per line, no header) to perform enrichment analysis.")

    # --- Step 1: File Upload (Stays at the top) ---
    uploaded_sample_file = st.file_uploader(
        "Upload your list of CpG sites", type=['csv', 'txt']
    )

    if uploaded_sample_file:
        # --- Step 2: Configure Analysis Parameters (New Layout) ---
        st.header("Configure Analysis Parameters")

        # Create two columns for the configuration options
        col1, col2 = st.columns(2)

        custom_background_file = None  # Initialize to None

        # --- Left Column: Background Selection ---
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
        
        # --- Right Column: Statistical Options ---
        with col2:
            st.subheader("Statistical Options")
            p_value_threshold = st.number_input(
                label="P-value threshold",
                min_value=0.0,
                max_value=1.0,
                value=0.05,  
                step=0.01,
                format="%.2f", 
                help="The significance threshold for the Fisher's Exact Test."
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
            
            # Save the new files for THIS run into the clean directory.
            user_sample_path = user_dir / USER_SAMPLE_NAME
            save_uploaded_file(uploaded_sample_file, user_sample_path)
            
            background = selected_background_name 
            if selected_background_name == 'custom':
                background = USER_CUSTOM_BACKGROUND_NAME 
                user_background_path = user_dir / USER_CUSTOM_BACKGROUND_NAME
                save_uploaded_file(custom_background_file, user_background_path)

            # Run the enrichment analysis pipeline
            run_enrichment_pipeline(user_dir, background, p_value_threshold)
            display_results()
            
    else:
        st.info("Please upload a file to begin.")


if __name__ == "__main__":
    main()