import streamlit as st
import pathlib
import shutil 
from statsmodels.stats.multitest import multitest_methods_names

from App_Utils.config import *
from App_Utils.utils import *
from App_Utils.pipeline import *
from App_Utils.view import *

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
                display_single_sample_results()
            
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

                # Call the merging function
                st.header("Comparative Analysis")
                comparison_results_dir = process_and_merge_comparison_results(user_dir)
                display_comparison_results(comparison_results_dir)

            
    else:
        st.info("Please provide the required sample file(s) to proceed.")


if __name__ == "__main__":
    main()  