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
        ("Single Sample Enrichment", "Two Sample Comparison", "Multi-Sample Group Comparison"),
        horizontal=True,
        label_visibility="collapsed"
    )
    
    # Initialize variables for sample sources
    sample_source_1 = None
    sample_source_2 = None
    sample_source_multi = None
    
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
            
            available_choices = get_background_options()
            selected_background_name_A = st.selectbox("Choose background for Sample A", options=available_choices, key="bg_A")

        with col_b:
            st.subheader("Sample B")
            uploaded_file_B = st.file_uploader("Upload file for Group B", type=['csv', 'txt'], key="uploader_B")
            if uploaded_file_B:
                sample_source_2 = uploaded_file_B

            selected_background_name_B = st.selectbox("Choose background for Sample B", options=available_choices, key="bg_B")

    elif analysis_type == "Multi-Sample Group Comparison":
        st.header("Multi-Sample Group Input")
        st.markdown("Upload a `.zip` file containing multiple sample files (one ID per line in each file).")
        uploaded_zip_file = st.file_uploader(
            "Upload your zip archive of CpG site lists", type=['zip'], key="multi_uploader"
        )
        if uploaded_zip_file:
            sample_source_multi = uploaded_zip_file
        
        group_comparison_method = st.selectbox(
            "Select Group Comparison Method",
            options=["fisher", "ttest"],
            index=1,  # Default to 'fisher'
            help="Choose the statistical method for comparing the real vs. control groups. 'fisher' compares enrichment frequency (recommended), while 'ttest' compares the average Odds Ratios."
        )

    # This part of the logic runs if at least one sample is provided 
    # The condition will need to be updated when we add the comparison logic
    if sample_source_1 or sample_source_multi: # For now, we only proceed if the first sample is there
        
        # Configure Analysis Parameters (This section remains mostly the same) 
        st.header("Configure Analysis Parameters")
        col1, col2, col3 = st.columns(3)
        selected_background_name = None
        custom_background_file = None
        
        # Make the background selector conditional
        if analysis_type == "Single Sample Enrichment" or analysis_type == "Multi-Sample Group Comparison":
            with col3:
                st.subheader("Background Universe")
                available_choices = get_background_options() + ['custom']
                selected_background_name = st.selectbox("Choose the background set", options=available_choices)
                if selected_background_name == 'custom':
                    st.info("Your custom background should be a file with one CpG ID per line.")
                    custom_background_file = st.file_uploader("Upload your custom background file", type=['csv', 'txt'], key='custom_bg')
        
        with col1:
            st.subheader("Statistical Options")
            p_value_threshold = st.number_input("P-value threshold", 0.0, 1.0, 0.05, 0.01, "%.2f")
        
        with col2:
            st.subheader("Correction Method")
            correction_options = sorted(list(multitest_methods_names.values()) + ['none'])
            default_index = correction_options.index('Bonferroni')
            selected_correction_method = st.selectbox("Multiple Testing Correction", options=correction_options, index=default_index)
        
        # Run Button Logic, check both samples for comparison mode
        run_button_disabled = (
            (analysis_type == "Single Sample Enrichment" and sample_source_1 is None) or
            (analysis_type == "Two Sample Comparison" and (sample_source_1 is None or sample_source_2 is None)) or
            (analysis_type == "Multi-Sample Group Comparison" and sample_source_multi is None) or
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
                display_single_sample_results(output_dir)
            
            elif analysis_type == "Two Sample Comparison":
                user_sample_a_path = user_dir / TOW_SAMPLE_COMPARISON_NAME_1
                user_sample_b_path = user_dir / TOW_SAMPLE_COMPARISON_NAME_2
                
                user_sample_a_path.mkdir(parents=True, exist_ok=True)
                user_sample_b_path.mkdir(parents=True, exist_ok=True)
                
                save_uploaded_file(sample_source_1, user_sample_a_path/ USER_SAMPLE_NAME)
                save_uploaded_file(sample_source_2, user_sample_b_path/ USER_SAMPLE_NAME)
                
                run_enrichment_pipeline_2samples(
                    user_dir=user_dir,
                    background_A=selected_background_name_A, 
                    background_B=selected_background_name_B, 
                    p_value=p_value_threshold,
                    correction_method=selected_correction_method
                )

                # Call the merging function
                st.header("Comparative Analysis")
                comparison_results_dir = process_and_merge_comparison_results(user_dir)
                display_comparison_results(comparison_results_dir)
            
            elif analysis_type == "Multi-Sample Group Comparison":
                # Define directories for real, control, and results
                real_samples_dir = user_dir / "real_samples"
                control_samples_dir = user_dir / "control_samples"
                multi_sample_results_dir = user_dir / "multi_sample_results"
                
                real_samples_dir.mkdir(parents=True, exist_ok=True)
                control_samples_dir.mkdir(parents=True, exist_ok=True)
                multi_sample_results_dir.mkdir(parents=True, exist_ok=True)

                # Unzip user-provided samples 
                handle_zip_extraction(sample_source_multi, real_samples_dir)
                organize_extracted_samples(real_samples_dir)

                # Generate control samples
                st.info("Generating control samples for comparison...")
                # Determine the full background path
                background_file_name = next((f.name for f in COMMON_BACKGROUND_ROOT.glob(f'{background}.*')), None)
                if background == "custom":
                    background_path = user_dir / USER_CUSTOM_BACKGROUND_NAME
                else:
                    background_path = COMMON_BACKGROUND_ROOT / background_file_name

                if not background_path.exists():
                    st.error(f"Could not find the specified background file: {background_path}")
                    st.stop()
                
                create_control_samples(real_samples_dir, background_path, control_samples_dir)

                # Run the parallel pipeline 
                run_enrichment_pipeline_multi_sample(
                    real_samples_dir=real_samples_dir,
                    control_samples_dir=control_samples_dir,
                    background=background,
                    p_value=p_value_threshold,
                    correction_method=selected_correction_method
                )

                # Set up the folder structure
                group_data_per_module(
                    real_samples_dir=real_samples_dir,
                    control_samples_dir=control_samples_dir,
                    multi_sample_results_dir=multi_sample_results_dir,
                )

                # Run the analysis in parallel
                run_modular_analysis_parallel(
                    multi_sample_results_dir=multi_sample_results_dir,
                    modules_dir=APP_DIR / "Modules",
                    method=group_comparison_method
                )
                
                display_multi_sample_results(multi_sample_results_dir)


            
    else:
        st.info("Please provide the required sample file(s) to proceed.")
    
    st.markdown("---")
    footer_html = """
    <div style="text-align: center; padding-top: 20px; color: grey;">
        <p>
            This is an open-source project. Check out the code on 
            <a href="https://github.com/Serotav/Enrichard" target="_blank" style="color: grey; text-decoration: underline;">
            GitHub
            </a>
        </p>
    </div>
    """
    st.markdown(footer_html, unsafe_allow_html=True)

if __name__ == "__main__":
    main()  