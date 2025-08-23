import streamlit as st
import subprocess
import pathlib
import polars as pl
import time
import shutil

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
            f"""Command: {' '.join(command)}\nReturn Code: {process.returncode}\n\n--- STDOUT ---\n{process.stdout}\n\n--- STDERR ---\n{process.stderr}""",
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
            f"""Command: {' '.join(command_a)}\nReturn Code: {returncode_a}\n\n--- STDOUT ---\n{stdout_a}\n\n--- STDERR ---\n{stderr_a}""",
            language='log'
        )
    
    with st.expander("Show Pipeline Execution Log for Sample B"):
        st.code(
            f"""Command: {' '.join(command_b)}\nReturn Code: {returncode_b}\n\n--- STDOUT ---\n{stdout_b}\n\n--- STDERR ---\n{stderr_b}""",
            language='log'
        )

    # Check for Success
    if returncode_a != 0 or returncode_b != 0:
        st.error("One or both parallel pipeline executions failed. Please check the logs above for details.")
        st.stop()
    else:
        st.success("Both pipelines completed successfully in parallel!")

def run_enrichment_pipeline_multi_sample(
    real_samples_dir: pathlib.Path, 
    control_samples_dir: pathlib.Path, 
    background: str, 
    p_value: float, 
    correction_method: str
) -> None:
    """
    Executes the enrichment pipeline for multiple real and control samples in parallel.
    """
    commands = []
    
    # Prepare commands for real samples
    for sample_path in real_samples_dir.iterdir():
        if not sample_path.is_dir(): continue
        output_dir = sample_path / USER_MODULE_OUTPUT
        output_dir.mkdir(parents=True, exist_ok=True)
        cmd = ['bash', str(MASTER_SCRIPT_PATH), str(sample_path), str(output_dir), str(background), str(p_value), str(correction_method)]
        commands.append({'cmd': cmd, 'name': f"Real: {sample_path.name}", 'log_expander': None})

    # Prepare commands for control samples
    for sample_path in control_samples_dir.iterdir():
        if not sample_path.is_dir(): continue
        output_dir = sample_path / USER_MODULE_OUTPUT
        output_dir.mkdir(parents=True, exist_ok=True)
        cmd = ['bash', str(MASTER_SCRIPT_PATH), str(sample_path), str(output_dir), str(background), str(p_value), str(correction_method)]
        commands.append({'cmd': cmd, 'name': f"Control: {sample_path.name}", 'log_expander': None})

    st.markdown("---")
    st.info(f"Starting parallel analysis for {len(commands)} samples...")

    # Create expanders for logs first


    with st.spinner(f"Running enrichment for {len(commands)} samples in parallel. This may take a while..."):
        processes = [(item, subprocess.Popen(item['cmd'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)) for item in commands]
        
        # Wait for all processes to complete
        for item, process in processes:
            stdout, stderr = process.communicate()
            item['stdout'] = stdout
            item['stderr'] = stderr
            item['returncode'] = process.returncode

    st.success("All parallel pipelines completed!")

    # Display logs and check for errors
    any_errors = False
    for item in commands:
        if item['returncode'] != 0:
            any_errors = True
            st.error(f"Pipeline for {item['name']} failed. Check log above.")

    if any_errors:
        st.error("One or more pipeline executions failed. Please check the logs above for details.")
        st.stop()


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

        if "_HEATMAP" in str(file_a):
            st.info(f"heatmap for {str(file_a)}")
            fake_heatmap = comparison_dir / module_name / f"{module_name}_HEATMAP.tsv"
            fake_heatmap.parent.mkdir(exist_ok=True)
            fake_heatmap.touch()
            continue
        
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

def group_data_per_module(
    real_samples_dir: pathlib.Path,
    control_samples_dir: pathlib.Path,
    multi_sample_results_dir: pathlib.Path,
) -> None:
    """
    Organizes results into a modular structure for group analysis.
    """
    st.markdown("---")
    st.info("Organizing results for group analysis...")

    # 1. Determine which modules were run
    try:
        first_sample_output = next(real_samples_dir.iterdir()) / USER_MODULE_OUTPUT
        module_names = [d.name for d in first_sample_output.iterdir() if d.is_dir()]
    except StopIteration:
        st.warning("No sample results found. Cannot proceed with group analysis.")
        return

    with st.spinner("Copying result files for group analysis..."):
        for module_name in module_names:
            # 2. Define paths for this module
            module_analysis_dir = multi_sample_results_dir / module_name
            module_real_results_dir = module_analysis_dir / "real_results"
            module_control_results_dir = module_analysis_dir / "control_results"
            
            module_analysis_dir.mkdir(parents=True, exist_ok=True)
            module_real_results_dir.mkdir(exist_ok=True)
            module_control_results_dir.mkdir(exist_ok=True)

            # 3. Gather and copy result files
            _gather_results_for_module(real_samples_dir, module_name, module_real_results_dir)
            _gather_results_for_module(control_samples_dir, module_name, module_control_results_dir)

    # We don't need those 2 folders anymore
    shutil.rmtree(real_samples_dir)
    shutil.rmtree(control_samples_dir)
    st.success("Successfully organized result files.")

def _gather_results_for_module(source_dir: pathlib.Path, module_name: str, target_dir: pathlib.Path):
    """
    Copies all result files for a specific module from the source sample directories to the target directory.
    """
    for sample_dir in source_dir.iterdir():
        if not sample_dir.is_dir(): continue
        
        # Path to the module's output directory for a single sample
        module_output_dir = sample_dir / USER_MODULE_OUTPUT / module_name
        
        if module_output_dir.exists():
            # Copy all .tsv files (handles both single file and heatmaps)
            for result_file in module_output_dir.glob('*.tsv'):
                shutil.copy(result_file, target_dir / f"{sample_dir.name}_{result_file.name}")

def run_modular_analysis_parallel(
    multi_sample_results_dir: pathlib.Path,
    modules_dir: pathlib.Path,
    method: str
) -> None:
    """
    Finds and executes the 'multisample.sh' script for each module in parallel.
    """
    st.markdown("---")
    st.info("Running final group analysis for each module in parallel...")

    commands = []
    module_dirs = [d for d in multi_sample_results_dir.iterdir() if d.is_dir()]

    for module_analysis_dir in module_dirs:
        module_name = module_analysis_dir.name
        analysis_script_path = modules_dir / module_name / "multisample.sh"
        
        if analysis_script_path.exists():
            command = ["bash", str(analysis_script_path), str(module_analysis_dir), method]
            commands.append({'cmd': command, 'name': module_name, 'log_expander': None})
        else:
            st.warning(f"Analysis script `multisample.sh` not found for module '{module_name}'. Skipping.")

    if not commands:
        st.warning("No analysis scripts found for any module.")
        return

    # Create expanders for logs first
    for item in commands:
        item['log_expander'] = st.expander(f"Show Group Analysis Log for {item['name']}", expanded=False)

    with st.spinner(f"Running group analysis for {len(commands)} modules..."):
        processes = [(item, subprocess.Popen(item['cmd'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)) for item in commands]
        
        for item, process in processes:
            stdout, stderr = process.communicate()
            item['stdout'] = stdout
            item['stderr'] = stderr
            item['returncode'] = process.returncode

    st.success("All group analyses completed!")

    # Display logs and check for errors
    any_errors = False
    for item in commands:
        with item['log_expander']:
            st.code(
                f"""Command: {' '.join(item['cmd'])}\nReturn Code: {item['returncode']}\n\n--- STDOUT ---\n{item['stdout']}\n\n--- STDERR ---\n{item['stderr']}""",
                language='log'
            )
        if item['returncode'] != 0:
            any_errors = True
            st.error(f"Group analysis for module '{item['name']}' failed. Check log above.")

    if any_errors:
        st.error("One or more group analyses failed.")
        st.stop()