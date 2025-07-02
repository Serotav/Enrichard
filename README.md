# Enrichard: CpG Site Enrichment Analysis Tool

**Enrichard** is a powerful and user-friendly web application designed for CpG site enrichment analysis. It allows researchers to identify over-represented biological traits, pathways, or genomic features within a given set of CpG probe IDs. The tool is built with Python, Streamlit, and various data analysis libraries, providing a seamless and interactive experience for both single-sample and two-sample comparative analyses.

The project is built with a **modular architecture**, making it easy to extend and customize. New analysis modules can be added with minimal effort, and existing ones can be modified to suit specific research needs. To ensure maximum performance and efficiency, the data processing pipeline is optimized with **Polars**, a lightning-fast data manipulation library.

## 🚀 Key Features

- **Single Sample Enrichment:** Analyze a single list of CpG sites to find enriched biological terms.
- **Two-Sample Comparison:** Compare two distinct sets of CpG sites to identify common and unique enriched traits between them.
- **Custom Backgrounds:** Use predefined genomic backgrounds or upload your own custom background for more tailored analyses.
- **Multiple Correction Methods:** Choose from a variety of statistical correction methods, including Bonferroni and FDR, to adjust p-values and reduce false positives.
- **Interactive Visualizations:** Explore your results through interactive tables and visualizations, making it easy to interpret complex data.
- **Parallel Processing:** Leverages parallel processing to deliver fast and efficient analysis, even with large datasets.
- **Open-Source:** The entire project is open-source, encouraging transparency, collaboration, and community-driven improvements.

## ⚙️ How It Works

Enrichard's backend is built on a robust and modular pipeline that performs the following steps:

1.  **Data Input:** Users can upload a file containing a list of CpG probe IDs or use one of the provided example datasets.
2.  **Enrichment Analysis:** The application uses Fisher's exact test to determine the statistical significance of the enrichment of each biological trait in the provided sample.
3.  **Multiple Testing Correction:** To account for the large number of tests performed, the raw p-values are adjusted using a selected multiple testing correction method.
4.  **Results Visualization:** The final results, including the enriched traits, p-values, and other relevant statistics, are presented in a clear and interactive table.

For two-sample comparisons, the pipeline runs the enrichment analysis for both samples in parallel and then identifies the common traits, providing a side-by-side comparison of the results.

## 🛠️ Setup and Installation

To run Enrichard locally, you will need to have Docker and Docker Compose installed.

1.  **Clone the repository:**

    ```bash
    git clone https://github.com/Serotav/Enrichard.git
    cd Enrichard
    ```

2.  **Build and run the Docker container:**

    ```bash
    docker-compose up --build
    ```

3.  **Access the application:**

    Open your web browser and navigate to `http://localhost:8501` to start using Enrichard.

