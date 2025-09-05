# Enrichard

**Enrichard**, also known as FLASH (Fast Loci Annotation of Significant enricHments), is a powerful and user-friendly web application designed for CpG site enrichment analysis. It allows researchers to identify over-represented biological traits, pathways, or genomic features within a given set of CpG probe IDs. The tool is built with Python, Streamlit, and a high-performance data processing backend, providing a seamless and interactive experience for complex epigenetic analyses.

The project is built with a **modular architecture**, making it easy to extend and customize. New analysis modules can be added with minimal effort. To ensure maximum performance and efficiency, the data processing pipeline is heavily optimized with **Polars**.

## 🚀 Features

- **Three Analisys Pipelines:**
    
    - **Single Sample Enrichment:** Analyze a single list of CpG sites to find enriched biological terms.
        
    - **Two-Sample Comparison:** Compare two distinct sets of CpG sites to identify common and unique enriched traits.
        
    - **Multi-Sample Group Comparison:** Perform a meta-analysis on entire cohorts to identify robust biological themes that are consistently enriched across many samples compared to random controls.

- **High Performance:** The analysis core is designed for speed, leveraging parallel processing, pre-computation, and the Polars library to deliver results in minutes, not hours.
    
- **Custom Backgrounds:** Use predefined genomic backgrounds (HM450, EPIC, EPIC+, EPICv2) or upload your own.
    
- **Multiple Correction Methods:** Choose from a variety of statistical correction methods, including Bonferroni and FDR (Benjamini-Hochberg), to control for false positives.
    
- **Interactive Visualizations:** Explore your results through dynamic dot plots, heatmaps, and dumbbell plots, making it easy to interpret complex data.

## ⚙️ How It Works

Enrichard's backend is built on a robust and modular pipeline that performs the following steps:

1.  **Data Input:** Users can upload a file containing a list of CpG probe IDs or use one of the provided example datasets.
2.  **Enrichment Analysis:** The application uses Fisher's exact test to determine the statistical significance of the enrichment of each biological trait in the provided sample.
3.  **Multiple Testing Correction:** To account for the large number of tests performed, the raw p-values are adjusted using a selected multiple testing correction method.
4.  **Results Visualization:** The final results, including the enriched traits, p-values, and other relevant statistics, are presented in a clear and interactive table.

For two-sample comparisons, the pipeline runs the enrichment analysis for both samples in parallel and then identifies the common traits, providing a side-by-side comparison of the results.

## 🛠️ Setup and Installation

To run Enrichard locally, you will need to have **Docker**, **Docker Compose**, and **Git LFS** installed.

Some modules, use large pre-computed data files that are stored using Git Large File Storage (LFS). 
```bash
# Install Git LFS
sudo pacman -S --needed git-lfs

# Install for your user account
git lfs install
```

```
### Clone the repository and ensure you pull the LFS files.
```bash
git clone https://github.com/Serotav/Enrichard
cd Enrichard
git lfs pull
docker compose up --build -d 
```