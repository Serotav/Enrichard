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

def run_enrichment_pipeline(user_dir: pathlib.Path, outputfoler:pathlib.Path, background: str, p_value: float, correction_method: str) -> None:
    """
    Executes the master enrichment shell script and displays its output.
    """
    command = ['bash', str(MASTER_SCRIPT_PATH), str(user_dir), str(outputfoler), str(background), str(p_value), str(correction_method)]
    
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

def run_enrichment_pipeline_2samples( user_dir: pathlib.Path, background: str, p_value: float,correction_method: str) -> None:
    """
    Executes the enrichment pipeline for two samples in PARALLEL.
    """
    # Define paths for sample A and sample B
    user_dir_module_A = user_dir / "Sample_A"
    user_dir_module_B = user_dir / "Sample_B"
    
    # Define output directories for each sample's results
    output_dir_a = user_dir_module_A / USER_MODULE_OUTPUT
    output_dir_b = user_dir_module_B / USER_MODULE_OUTPUT
    output_dir_a.mkdir(parents=True, exist_ok=True)
    output_dir_b.mkdir(parents=True, exist_ok=True)

    # Commands for both pipelines
    command_a = [
        'bash', str(MASTER_SCRIPT_PATH),
        str(user_dir_module_A), str(output_dir_a),
        str(background), str(p_value), str(correction_method)
    ]
    command_b = [
        'bash', str(MASTER_SCRIPT_PATH),
        str(user_dir_module_B), str(output_dir_b),
        str(background), str(p_value), str(correction_method)
    ]

    st.markdown("---")
    st.info("Starting parallel analysis for two samples...")

    with st.spinner("Running enrichment for Sample A and Sample B in parallel. We'll be right back with results!"):
        # Launch both processes using Popen without waiting 
        # stdout and stderr are piped so we can capture them later.
        process_a = subprocess.Popen(command_a, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        process_b = subprocess.Popen(command_b, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        # Now, wait for both processes to complete and capture their output
        # The .communicate() method waits for the process to finish and reads all output.
        stdout_a, stderr_a = process_a.communicate()
        stdout_b, stderr_b = process_b.communicate()

        # We get the return codes after they have finished.
        returncode_a = process_a.returncode
        returncode_b = process_b.returncode

    # Display Logs (logic is now outside the spinner) ---
    with st.expander("Show Pipeline Execution Log for Sample A"):
        st.code(
            f"Command: {' '.join(command_a)}\nReturn Code: {returncode_a}\n\n"
            f"--- STDOUT ---\n{stdout_a}\n\n--- STDERR ---\n{stderr_a}",
            language='log'
        )
    
    with st.expander("Show Pipeline Execution Log for Sample B"):
        st.code(
            f"Command: {' '.join(command_b)}\nReturn Code: {returncode_b}\n\n"
            f"--- STDOUT ---\n{stdout_b}\n\n--- STDERR ---\n{stderr_b}",
            language='log'
        )

    # --- Check for Success ---
    if returncode_a != 0 or returncode_b != 0:
        st.error("One or both parallel pipeline executions failed. Please check the logs above for details.")
        st.stop()
    else:
        st.success("Both pipelines completed successfully in parallel!")

def process_and_merge_comparison_results(output_base_dir: pathlib.Path) -> pathlib.Path:
    """
    Finds common results between Sample A and Sample B for each module.
    Returns: pathlib.Path: The path to the directory containing the merged results.
    """
    dir_a_modules = output_base_dir / "Sample_A" / USER_MODULE_OUTPUT
    dir_b_modules = output_base_dir / "Sample_B" / USER_MODULE_OUTPUT
    comparison_dir = output_base_dir / "Comparison_Results"
    comparison_dir.mkdir(exist_ok=True)

    st.info("Processing and merging results from Sample A and Sample B...")
    
    # Find all result files in Sample A's directory structure
    result_files_a = list(dir_a_modules.glob('*/*.tsv'))
    
    if not result_files_a:
        st.warning(f"No result files found for Sample A {dir_a_modules} {result_files_a}. Cannot perform comparison.")
        return None

    for file_a in result_files_a:
        file_b = dir_b_modules / file_a.parent.name / file_a.name
        module_name = file_a.parent.name
        
        if not file_b.exists():
            st.warning(f"Result file for module '{module_name}' not found for Sample B. Skipping.")
            continue
            
        try:
            df_a = pl.read_csv(file_a, separator='\t')
            df_b = pl.read_csv(file_b, separator='\t')
            
            # Perform the Inner Join to find common traits 
            # We add suffixes to distinguish columns from A and B after the join
            merged_df = df_a.join(
                df_b,
                on="Trait",
                how="inner",
                suffix="_B"
            )

            # Save the merged result 
            output_module_dir = comparison_dir / module_name
            output_module_dir.mkdir(exist_ok=True)
            output_file_path = output_module_dir / f"{module_name}_common.tsv"
            
            merged_df.write_csv(output_file_path, separator='\t')

        except Exception as e:
            st.error(f"Failed to process and merge results for module '{module_name}': {e}")
    
    st.success("Comparison processing complete!")
    return comparison_dir

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

def create_dumbbell_plot(df: pd.DataFrame):
    """
    Creates an Altair dumbbell plot to compare enrichment results from two samples.
    """
    # Keep top 20 traits, sorted by the most significant p-value between the two samples
    df_to_plot = df.assign(
        min_P_Value=df[['P-Value', 'P-Value_B']].min(axis=1)
    ).sort_values('min_P_Value').head(20).copy()

    sort_order = df_to_plot['Trait'].tolist()

    # Base chart for common encodings
    base = alt.Chart(df_to_plot).encode(
        y=alt.Y('Trait:N', sort=sort_order, title="Enriched Trait"),
        tooltip=[
            alt.Tooltip('Trait:N'),
            alt.Tooltip('Fold-Change:Q', title="Fold Change (A)", format=".2f"),
            alt.Tooltip('Fold-Change_B:Q', title="Fold Change (B)", format=".2f"),
            alt.Tooltip('P-Value:Q', title="P-Value (A)", format=".2e"),
            alt.Tooltip('P-Value_B:Q', title="P-Value (B)", format=".2e"),
        ]
    )

    # The connecting line (the "bar" of the dumbbell)
    line = base.mark_rule().encode(
        x=alt.X('Fold-Change:Q', title="Fold Change", scale=alt.Scale(zero=False)),
        x2=alt.X2('Fold-Change_B:Q'),
    )

    # The dots for Sample A
    points_a = base.mark_circle(size=100, color='#1f77b4').encode( # Blue
        x=alt.X('Fold-Change:Q'),
        size=alt.Size('a:Q', legend=alt.Legend(title="Count in Sample"))
    )
    
    # The dots for Sample B
    points_b = base.mark_circle(size=100, color='#ff7f0e').encode( # Orange
        x=alt.X('Fold-Change_B:Q'),
        size=alt.Size('a_B:Q') # Legend is shared with points_a
    )

    # Layer the three charts together
    chart = (line + points_a + points_b).properties(
        title="Comparison of Common Enriched Traits"
    ).interactive()
    
    return chart

def display_comparison_results(comparison_dir: pathlib.Path):
    """
    Scans for merged result files and displays them in a comparative view.
    """
    result_files = sorted(list(comparison_dir.glob('*/*.tsv')))
    
    if not result_files:
        st.warning("Comparison processing complete, but no common enriched traits were found in any module.")
        return

    st.header("Comparative Enrichment Results")
    st.markdown("This view shows traits that were significantly enriched in **both** Sample A and Sample B.")
    st.markdown("---")

    module_names = [path.parent.name for path in result_files]
    tabs = st.tabs(module_names)

    # Define paths to the original single-sample results
    output_base_dir = comparison_dir.parent
    dir_a = output_base_dir / "Sample_A" / USER_MODULE_OUTPUT
    dir_b = output_base_dir / "Sample_B" / USER_MODULE_OUTPUT

    for i, file_path in enumerate(result_files):
        with tabs[i]:
            module_name = module_names[i]
            st.markdown(f"#### Module: `{module_name}`")
            
            try:
                # Load the merged (common) results
                merged_df = pl.read_csv(file_path, separator='\t').to_pandas()
                
                if merged_df.empty:
                    st.info("No common significant traits for this module in merge.")

                # --- Create and display the comparison chart ---
                st.subheader("Comparison Plot")
                comparison_chart = create_dumbbell_plot(merged_df)
                st.altair_chart(comparison_chart, use_container_width=True)

                # --- Display the data tables in expanders ---
                st.subheader("Data Tables")
                
                with st.expander("Show Common Results Data (Merged Table)"):
                    # Rename columns for clarity before displaying
                    display_df = merged_df.rename(columns={
                        'P-Value': 'P-Value_A', 'Fold-Change': 'Fold-Change_A', 'a': 'Count_A',
                        'b': 'b_A', 'c': 'c_A', 'd': 'd_A',
                        'a_B': 'Count_B'
                    })
                    float_cols = display_df.select_dtypes(include='float').columns
                    format_dict = {col: '{:.2e}' for col in float_cols}
                    st.dataframe(display_df.style.format(format_dict))

                # Find and load the original single-sample data files
                original_file_a = dir_a / module_name / f"{module_name}.tsv"
                original_file_b = dir_b / module_name / f"{module_name}.tsv"

                with st.expander("Show Full Results for Sample A"):
                    if original_file_a.exists():
                        df = pl.read_csv(original_file_a, separator='\t').to_pandas()
                        float_cols = df.select_dtypes(include='float').columns
                        format_dict = {col: '{:.2e}' for col in float_cols}
                        st.dataframe(df.style.format(format_dict))
                    else:
                        st.warning("Original result file for Sample A not found.")
                
                with st.expander("Show Full Results for Sample B"):
                    if original_file_b.exists():
                        df = pl.read_csv(original_file_b, separator='\t').to_pandas()
                        float_cols = df.select_dtypes(include='float').columns
                        format_dict = {col: '{:.2e}' for col in float_cols}
                        st.dataframe(df.style.format(format_dict))
                    else:
        
                        st.warning("Original result file for Sample B not found.")
            
            except Exception as e:
                st.error(f"Error displaying comparison results for module '{module_name}': {e}")


def main():
    st.set_page_config(layout="wide")

    st.markdown("<h1 style='text-align: center;'>🚀 E N R I C H A R D 🚀</h1>", unsafe_allow_html=True)
    st.subheader("CpG Site Enrichment Analysis")
    initialize_session()

    st.markdown("---")

    analysis_type = st.radio(
        "Select Analysis Type",
        ("Single Sample Enrichment", "Two Sample Comparison"),
        horizontal=True,
        label_visibility="collapsed"
    )
    
    # Initialize variables for sample sources
    sample_source_1 = None
    sample_source_2 = None
    
    # Section for Single Sample Analysis
    if analysis_type == "Single Sample Enrichment":
        st.header("Single Sample Input")
        st.markdown("Upload a file containing CpG probe IDs (one ID per line) to perform enrichment analysis.")
        
        input_method = st.radio(
            "Choose your input method:",
            ("Upload a File", "Use an Example"),
            horizontal=True,
            label_visibility="collapsed",
            key="single_sample_input_method" # Use a key to avoid widget conflicts
        )

        if input_method == "Upload a File":
            uploaded_sample_file = st.file_uploader(
                "Upload your list of CpG sites", type=['csv', 'txt'], key="single_uploader"
            )
            if uploaded_sample_file:
                sample_source_1 = uploaded_sample_file

        elif input_method == "Use an Example":
            example_files = get_example_files()
            if not example_files:
                st.info("No example files found.")
            else:
                options = ["- Select an example -"] + list(example_files.keys())
                selected_example_name = st.selectbox(
                    "Choose an example file to run:",
                    options=options,
                    key="single_example"
                )
                if selected_example_name != "- Select an example -":
                    sample_source_1 = example_files[selected_example_name]

    #  Section for Two Sample Comparison 
    elif analysis_type == "Two Sample Comparison":
        st.header("Two Sample Input")
        st.markdown("Upload two files, each containing a list of CpG probe IDs. The analysis will identify common enriched traits.")

        col_a, col_b = st.columns(2)
        
        with col_a:
            st.subheader("Sample A")
            uploaded_file_A = st.file_uploader("Upload file for Group A", type=['csv', 'txt'], key="uploader_A")
            if uploaded_file_A:
                sample_source_1 = uploaded_file_A

        with col_b:
            st.subheader("Sample B")
            uploaded_file_B = st.file_uploader("Upload file for Group B", type=['csv', 'txt'], key="uploader_B")
            if uploaded_file_B:
                sample_source_2 = uploaded_file_B

    # This part of the logic runs if at least one sample is provided 
    # The condition will need to be updated when we add the comparison logic
    if sample_source_1: # For now, we only proceed if the first sample is there
        
        # Configure Analysis Parameters (This section remains mostly the same) 
        st.header("Configure Analysis Parameters")
        col1, col2, col3 = st.columns(3)
        custom_background_file = None
        
        with col1:
            st.subheader("Background Universe")
            available_choices = get_background_options() + ['custom']
            selected_background_name = st.selectbox("Choose the background set", options=available_choices)
            if selected_background_name == 'custom':
                st.info("Your custom background should be a file with one CpG ID per line.")
                custom_background_file = st.file_uploader("Upload your custom background file", type=['csv', 'txt'], key='custom_bg')
        
        with col2:
            st.subheader("Statistical Options")
            p_value_threshold = st.number_input("P-value threshold", 0.0, 1.0, 0.05, 0.01, "%.2f")
        
        with col3:
            st.subheader("Correction Method")
            correction_options = sorted(list(multitest_methods_names.values()) + ['none'])
            default_index = correction_options.index('Bonferroni')
            selected_correction_method = st.selectbox("Multiple Testing Correction", options=correction_options, index=default_index)
        
        # Run Button Logic, check both samples for comparison mode
        run_button_disabled = (
            (analysis_type == "Single Sample Enrichment" and sample_source_1 is None) or
            (analysis_type == "Two Sample Comparison" and (sample_source_1 is None or sample_source_2 is None)) or
            (selected_background_name == 'custom' and custom_background_file is None)
        )

        st.markdown("---") 

        if st.button("Run Analysis", type="primary", disabled=run_button_disabled):
            # Clean up and run the existing single-sample pipeline
            user_dir = st.session_state.user_dir
            if user_dir.exists(): shutil.rmtree(user_dir)
            user_dir.mkdir(exist_ok=True)
            
            user_sample_path = user_dir / USER_SAMPLE_NAME 
            output_dir = user_dir / USER_MODULE_OUTPUT

            # custom background handling
            background = selected_background_name
            if selected_background_name == 'custom':
                background = USER_CUSTOM_BACKGROUND_NAME
                user_background_path = user_dir / USER_CUSTOM_BACKGROUND_NAME
                save_uploaded_file(custom_background_file, user_background_path)

            if analysis_type == "Single Sample Enrichment":
                # Logic to handle both uploaded and example files for sample 1
                if isinstance(sample_source_1, pathlib.Path): # It's an example
                    st.info(f"Using example file: {sample_source_1.name}")
                    copy_example_file(sample_source_1, user_sample_path)
                else: # It's an uploaded file
                    st.info(f"Using uploaded file: {sample_source_1.name}")
                    save_uploaded_file(sample_source_1, user_sample_path)

                background = selected_background_name
                if selected_background_name == 'custom':
                    background = USER_CUSTOM_BACKGROUND_NAME
                    save_uploaded_file(custom_background_file, user_dir / USER_CUSTOM_BACKGROUND_NAME)

                run_enrichment_pipeline(user_dir, output_dir, background, p_value_threshold, selected_correction_method)
                display_results()
            
            elif analysis_type == "Two Sample Comparison":
                user_sample_a_path = user_dir / "Sample_A"
                user_sample_b_path = user_dir / "Sample_B"
                
                user_sample_a_path.mkdir(parents=True, exist_ok=True)
                user_sample_b_path.mkdir(parents=True, exist_ok=True)
                
                save_uploaded_file(sample_source_1, user_sample_a_path/ USER_SAMPLE_NAME)
                save_uploaded_file(sample_source_2, user_sample_b_path/ USER_SAMPLE_NAME)
                
                # Run the new two-sample pipeline function
                run_enrichment_pipeline_2samples(
                    user_dir=user_dir,
                    background=background,
                    p_value=p_value_threshold,
                    correction_method=selected_correction_method
                )

                # --- NEW: Call the merging function ---
                st.header("Comparative Analysis")
                comparison_results_dir = process_and_merge_comparison_results(user_dir)
                display_comparison_results(comparison_results_dir)

            
    else:
        st.info("Please provide the required sample file(s) to proceed.")


if __name__ == "__main__":
    main()