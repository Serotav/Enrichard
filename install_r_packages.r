# Install BiocManager if it's not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "http://cran.us.r-project.org")
}

# List of CRAN packages
cran_packages <- c("tidyverse")

# List of Bioconductor packages
bioc_packages <- c("AnnotationDbi", "org.Hs.eg.db")

# Install CRAN packages
install.packages(cran_packages, repos = "http://cran.us.r-project.org")

# Install Bioconductor packages
BiocManager::install(bioc_packages)