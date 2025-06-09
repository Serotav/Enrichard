if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "http://cran.us.r-project.org")
}

install.packages("tidyverse")

bioc_packages <- c("AnnotationDbi", "org.Hs.eg.db")
BiocManager::install(bioc_packages)