import streamlit as st
import subprocess
import pathlib
import polars as pl

from .config import *

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
    user_dir_module_A = user_dir / TOW_SAMPLE_COMPARISON_NAME_1
    user_dir_module_B = user_dir / TOW_SAMPLE_COMPARISON_NAME_2
    
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

    # Check for Success
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
    dir_a_modules = output_base_dir / TOW_SAMPLE_COMPARISON_NAME_1 / USER_MODULE_OUTPUT
    dir_b_modules = output_base_dir / TOW_SAMPLE_COMPARISON_NAME_2 / USER_MODULE_OUTPUT
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