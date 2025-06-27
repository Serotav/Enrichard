import streamlit as st
import pandas as pd
import altair as alt
import numpy as np
import polars as pl
import pathlib

from .config import *

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

    alt.Chart(df_to_plot).mark_circle().encode(
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
    
def display_single_sample_results()->None:
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
