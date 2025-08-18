import streamlit as st
import pandas as pd
import altair as alt
import numpy as np
import polars as pl
import pathlib

from .config import *

def display_altair_heatmap(file_path: pathlib.Path):
    """
    Reads a TSV file and displays it as an Altair heatmap with a dynamic color scale.
    """
    try:
        df = pl.read_csv(file_path, separator='\t')
        
        # Altair requires data in long format, so we melt the dataframe
        long_df = df.melt(
            id_vars=df.columns[0], 
            variable_name="State", 
            value_name="Odds-Ratio"
        ).rename({df.columns[0]: "Cell_Type"})

        # Determine the dynamic range for the color scale
        min_or = long_df["Odds-Ratio"].min()
        max_or = long_df["Odds-Ratio"].max()

        # Create a robust diverging domain centered at 1.0
        color_domain = [min_or, 1.0, max_or]
        if min_or == 1.0 and max_or == 1.0:
            color_domain = [0.9, 1.0, 1.1] # Handle edge case where all values are 1.0

        # Get the order of states for proper sorting on the x-axis
        state_order = df.columns[1:]
        cell_size, chart_height, chart_width= 15, 600, 800
        heatmap = alt.Chart(long_df).mark_rect(size=cell_size).encode( 
            x=alt.X('State:N', sort=state_order, title="Chromatin State"),
            y=alt.Y('Cell_Type:N', title="Cell Type", axis=alt.Axis(labelLimit=200)),
            color=alt.Color('Odds-Ratio:Q',
                scale=alt.Scale(scheme='redblue', domain=color_domain, reverse=True),
                legend=alt.Legend(title="Odds Ratio")
            ),
            tooltip=[
                alt.Tooltip('Cell_Type:N', title="Cell Type"),
                alt.Tooltip('State:N', title="State"),
                alt.Tooltip('Odds-Ratio:Q', title="Odds Ratio", format=".3f")
            ]
        ).properties(
            title="Chromatin State Enrichment Heatmap",
        ).interactive()

        st.altair_chart(heatmap, use_container_width=True)

    except Exception as e:
        st.error(f"Failed to create heatmap: {e}")

def create_dot_plot(df: pd.DataFrame)->None:
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

    return alt.Chart(df_to_plot).mark_circle().encode(
        y=alt.Y('Trait:N', sort=sort_order, title="Enriched Trait", axis=alt.Axis(labelLimit=200)),
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

def create_dumbbell_plot(df: pd.DataFrame)->None:
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


def render_methodology_explanation()->None:
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
        # Check if the file is a heatmap
        if "_HEATMAP" in file_path.name:
            st.subheader("Visualization")
            display_altair_heatmap(file_path)
            with st.expander("Show Full Data Table"):
                df = pl.read_csv(file_path, separator='\t').to_pandas()
                st.dataframe(df)
            return True

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
        st.error(f"Error reading or processing result file '{file_path.name} {file_path}': {e}")
        return False
    
def display_single_sample_results(module_output:str)->None:
    """
    Scans for result files, creates tabs, and displays results and methodology.
    """
    module_output = pathlib.Path(module_output)
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

def display_comparison_results(comparison_dir: pathlib.Path)->None:
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
    dir_a = output_base_dir / TOW_SAMPLE_COMPARISON_NAME_1 / USER_MODULE_OUTPUT
    dir_b = output_base_dir / TOW_SAMPLE_COMPARISON_NAME_2 / USER_MODULE_OUTPUT

    for i, file_path in enumerate(result_files):
        with tabs[i]:
            module_name = module_names[i]
            st.markdown(f"#### Module: `{module_name}`")
            
            try:
                # Handle heatmap files differently
                if "_HEATMAP" in file_path.name:
                    st.subheader("Heatmap Comparison")
                    st.info("Heatmap comparison view is not yet implemented. Showing individual heatmaps.")
                    
                    col1, col2 = st.columns(2)
                    with col1:
                        st.markdown("##### Sample A")
                        original_file_a = dir_a / module_name / file_path.name
                        if original_file_a.exists():
                            display_altair_heatmap(original_file_a)
                        else:
                            st.warning("Result file for Sample A not found.")
                    
                    with col2:
                        st.markdown("##### Sample B")
                        original_file_b = dir_b / module_name / file_path.name
                        if original_file_b.exists():
                            display_altair_heatmap(original_file_b)
                        else:
                            st.warning("Result file for Sample B not found.")
                    continue

                # Load the merged (common) results for standard dot plots
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

# --- MULTI SAMPLE ---
def create_group_dumbbell_plot(df: pd.DataFrame) -> alt.Chart:
    """
    Creates an Altair dumbbell plot for multi-sample group comparison results.
    Returns: An Altair Chart object representing the dumbbell plot.
    """
    if df.empty:
        st.info("The result data is empty; cannot generate a plot.")
        return None

    # --- 1. Detect analysis method to set dynamic chart elements ---
    method = df['Statistic_Name'].iloc[0]
    if method == "T-statistic":
        x_axis_title_base = "Mean log(Odds-Ratio)"
        statistic_tooltip_title = "T-statistic"
        chart_title = "Group Enrichment Comparison (T-test)"
    elif method == "Group_Odds_Ratio":
        x_axis_title_base = "Proportion of Enriched Samples"
        statistic_tooltip_title = "Group Odds Ratio"
        chart_title = "Group Enrichment Comparison (Fisher's Test)"
    else:
        x_axis_title_base = "Value"
        statistic_tooltip_title = "Statistic"
        chart_title = "Group Enrichment Comparison"

    # --- 2. Filter data and prepare for plotting ---
    df_to_plot = df.sort_values("P-value").head(20).copy()
    sort_order = df_to_plot['Trait'].tolist()

    # --- 3. Create the dumbbell plot components ---
    
    # The connecting line (the "bar") uses the original wide-format data
    line = alt.Chart(df_to_plot).mark_rule().encode(
        y=alt.Y('Trait:N', sort=sort_order, title="Enriched Trait"),
        x=alt.X('Value_Control:Q', title=f"{x_axis_title_base} (Value per Group)"),
        x2=alt.X2('Value_Real:Q'),
        tooltip=[
            alt.Tooltip('Trait:N', title="Trait"),
            alt.Tooltip('P-value:Q', title="P-Value", format=".2e"),
            alt.Tooltip('Statistic:Q', title=statistic_tooltip_title, format=".3f"),
        ]
    )

    # For the points and legend, we transform the data to a long format
    df_long = df_to_plot.melt(
        id_vars=['Trait', 'P-value', 'Statistic'],
        value_vars=['Value_Real', 'Value_Control'],
        var_name='Group',
        value_name='Value'
    ).replace({'Value_Real': 'Real', 'Value_Control': 'Control'})

    # The dots for both groups, colored by the new 'Group' column to create a legend
    points = alt.Chart(df_long).mark_circle(size=100, opacity=0.9).encode(
        y=alt.Y('Trait:N', sort=sort_order),
        x=alt.X('Value:Q'),
        color=alt.Color('Group:N',
            scale=alt.Scale(
                domain=['Real', 'Control'],
                range=['#1f77b4', '#ff7f0e'] # Blue for Real, Orange for Control
            ),
            legend=alt.Legend(title="Sample Group")
        ),
        tooltip=[
            alt.Tooltip('Trait:N', title="Trait"),
            alt.Tooltip('Group:N', title="Group"),
            alt.Tooltip('Value:Q', title="Value", format=".3f"),
            alt.Tooltip('P-value:Q', title="P-Value", format=".2e"),
            alt.Tooltip('Statistic:Q', title=statistic_tooltip_title, format=".3f"),
        ]
    )

    # --- 4. Layer the charts and add final properties ---
    chart = (line + points).properties(
        title=chart_title
    ).interactive()

    return chart

def display_multi_sample_results(multi_sample_results_dir: pathlib.Path):
    """
    Scans for group comparison result files, creates tabs for each module,
    and displays the results using the appropriate chart based on the filename.
    This function acts as a router, calling different plotting functions based on module type.
    """
    # General glob pattern to find all group comparison files
    result_files = sorted(list(multi_sample_results_dir.glob('*/group_comparison_*.tsv')))
    
    if not result_files:
        st.warning("Multi-sample group analysis complete, but no result files were found.")
        return

    st.header("Multi-Sample Group Comparison Results")

    with st.expander("How to Interpret These Results (Methodology)", expanded=False):
        render_group_analysis_methodology()
    
    st.markdown("---") 

    module_names = sorted(list(set([path.parent.name for path in result_files])))
    tabs = st.tabs(module_names)
    results_map = {path.parent.name: path for path in result_files}

    for i, module_name in enumerate(module_names):
        with tabs[i]:
            file_path = results_map[module_name]
            st.markdown(f"#### Module: `{module_name}`")
            
            try:
                df = pl.read_csv(file_path, separator='\t').to_pandas()
                if df.empty:
                    st.info("No significant group-level traits were found for this module.")
                    continue

                # --- ROUTER LOGIC: Choose plot based on filename ---
                
                # Case 1: Standard modules using dumbbell plot
                if "_dumbbell_" in file_path.name:
                    st.subheader("Group Enrichment Comparison")
                    chart = create_group_dumbbell_plot(df)
                    if chart:
                        st.altair_chart(chart, use_container_width=True)
                
                # Case 2: Chromatin module using heatmap plot
                elif "_heatmap_" in file_path.name:
                    st.subheader("Group Enrichment Heatmap")
                    chart = display_group_comparison_heatmap(df)
                    if chart:
                        st.altair_chart(chart, use_container_width=True)

                else:
                    st.warning(f"Could not determine plot type for '{file_path.name}'. Displaying raw data.")

                # Always show the data table in an expander
                with st.expander("Show Full Group Comparison Data"):
                    float_cols = df.select_dtypes(include='float').columns
                    format_dict = {col: '{:.2e}' for col in float_cols}
                    st.dataframe(df.style.format(format_dict))

            except Exception as e:
                st.error(f"Error displaying results for module '{module_name}': {e}")

def display_group_comparison_heatmap(df: pd.DataFrame) -> alt.Chart:
    """
    Creates an Altair heatmap from group-level chromatin analysis results.

    The color of each cell represents the T-statistic from the group comparison,
    indicating both the significance and direction of the difference between the
    real and control groups.

    Args:
        df: A pandas DataFrame containing the long-format group comparison results.
            Expected columns: 'Cell_Type', 'State', 'P-value', 'T-statistic', etc.

    Returns:
        An Altair Chart object representing the summary heatmap.
    """
    if df.empty:
        st.info("The result data is empty; cannot generate a plot.")
        return None

    # --- 1. Determine the dynamic range for the diverging color scale ---
    # We use the T-statistic: positive means Real > Control, negative means Control > Real.
    max_abs_t_stat = df['T-statistic'].abs().max()
    
    # Create a symmetrical domain around 0 for a balanced color scale
    color_domain = [-max_abs_t_stat, 0, max_abs_t_stat]
    
    # Handle the edge case where all statistics are zero
    if max_abs_t_stat == 0:
        color_domain = [-1, 0, 1]

    # --- 2. Build the Altair Heatmap ---
    heatmap = alt.Chart(df).mark_rect().encode(
        x=alt.X('State:N', title="Chromatin State", sort=None), # Use 'sort=None' to respect original order if possible
        y=alt.Y('Cell_Type:N', title="Cell Type", sort=alt.Sort(field="P-value", op="min")), # Sort rows by most significant P-value
        
        # Color is based on the T-statistic for direction and magnitude
        color=alt.Color('T-statistic:Q',
            scale=alt.Scale(scheme='redblue', domain=color_domain, reverse=True),
            legend=alt.Legend(title="T-statistic (Real vs Control)")
        ),
        
        # Tooltip provides the full details for each cell
        tooltip=[
            alt.Tooltip('Cell_Type:N', title="Cell Type"),
            alt.Tooltip('State:N', title="State"),
            alt.Tooltip('P-value:Q', title="Group P-Value", format=".2e"),
            alt.Tooltip('T-statistic:Q', title="Group T-statistic", format=".3f"),
            alt.Tooltip('Mean_logOR_Real:Q', title="Mean log(OR) Real", format=".3f"),
            alt.Tooltip('Mean_logOR_Control:Q', title="Mean log(OR) Control", format=".3f"),
        ]
    ).properties(
        title="Group Comparison of Chromatin State Enrichment"
    ).interactive()

    return heatmap


def render_group_analysis_methodology():
    """
    Renders a detailed explanation of the multi-sample group comparison methodologies,
    reflecting that the P-value cutoff is applied *after* the group-level test.
    """
    st.markdown("""
    The goal of this analysis is to identify traits that are consistently and significantly enriched across a group of your **Real Samples** when compared to a group of randomly generated **Control Samples**. This helps distinguish true biological signals from random statistical noise. Two different statistical methods are provided to perform this comparison.
    """)

    # Use tabs to cleanly separate the explanation for each method
    tab1, tab2 = st.tabs(["Fisher's Method (Frequency)", "T-test Method (Magnitude)"])

    with tab1:
        st.header("Fisher's Method: Comparing Frequencies")
        st.markdown("""
        This method answers the question: **"Is a trait found more *frequently* in the real samples than in the control samples?"**
        
        #### How It Works:
        1.  For every unique trait, the analysis counts how many of your **Real Sample** files contain a result for that trait. This count becomes `N_real_enriched`.
        2.  It does the same for the **Control Samples**, yielding `N_control_enriched`.
        3.  These counts are used to build a 2x2 contingency table for a group-level comparison:
        """)

        st.markdown("""
        | | Trait Present | Trait Absent |
        |---|---|---|
        | **Real Samples** | `N_real_enriched` | `N_real_total - N_real_enriched` |
        | **Control Samples**| `N_control_enriched`| `N_control_total - N_control_enriched` |
        """)
        
        st.markdown("""
        4. A **Fisher's Exact Test** is performed on this table to generate a new, **group-level P-value**.
        5. **Filtering:** Finally, the results from all tested traits are filtered. Only traits where this new **group-level P-value** is less than your chosen cutoff are shown in the final output.
        
        #### How to Interpret the Plot:
        - **X-Axis (`Proportion of Enriched Samples`):** Shows the fraction of samples in a group where the trait was present in the results.
        - **Blue Dot:** Represents the proportion for the **Real** sample group.
        - **Orange Dot:** Represents the proportion for the **Control** sample group.
        - **What to look for:** You are looking for traits where the **blue dot is significantly to the right** of the orange dot. This indicates a higher frequency of occurrence in your experimental samples.
        """)

    with tab2:
        st.header("T-test Method: Comparing Magnitudes")
        st.markdown("""
        This method answers the question: **"Is the *average strength* of enrichment (measured by log Odds-Ratio) significantly higher in the real group than in the control group?"**

        #### How It Works:
        1.  For each trait, the analysis collects the calculated **Odds Ratios (ORs)** from every real and control sample where it appeared.
        2.  **Imputation:** To create complete datasets, if a trait was not found in a sample's result file, we assign it an Odds Ratio of **1.0** (which corresponds to a log(OR) of **0**), representing a baseline of "no enrichment."
        3.  A **Welch's t-test** is performed on the two complete sets of log(Odds-Ratios) to generate a new, **group-level P-value**.
        4.  **Filtering:** Finally, the results from all tested traits are filtered. Only traits where this new **group-level P-value** from the t-test is less than your chosen cutoff are shown in the final output.

        #### How to Interpret the Plot:
        - **X-Axis (`Mean log(Odds-Ratio)`):** Shows the average enrichment strength for each group. A higher value means stronger enrichment on average.
        - **Blue Dot:** Represents the mean log(OR) for the **Real** sample group.
        - **Orange Dot:** Represents the mean log(OR) for the **Control** sample group.
        - **What to look for:** You are looking for traits where the **blue dot is significantly to the right** of the orange dot, indicating a higher average enrichment magnitude in your experimental samples.
        """)