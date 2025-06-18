if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "http://cran.us.r-project.org")
}

install.packages("tidyverse")

bioc_packages <- c("GenomicRanges", "TxDb.Hsapiens.UCSC.hg38.knownGene")
BiocManager::install(bioc_packages)