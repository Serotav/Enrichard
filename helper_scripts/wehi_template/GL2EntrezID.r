# This script reads a file with genomic coordinates and maps each location
# to its corresponding Entrez Gene ID by finding overlapping genes.
# Result: Probe_ID         EntrezID
# 0       cg00000029       5934
# Each probe_id might map to multiple Entrez IDs, so they might be duplicated.

# R needs to SHUT THE FUCK UP
library(tidyverse)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene) # Using hg38 version

# --- Configuration ---
args <- commandArgs(trailingOnly = TRUE)
INPUT_COORDINATES_FILE <- args[1]
OUTPUT_FILE <- args[2]
# Suppress package loading messages and warnings
# --- Step 1: Load Input Data and Create GRanges Object ---
message(paste("Loading coordinates from:", INPUT_COORDINATES_FILE))

# Read the input file and select the necessary columns
input_data <- read_tsv(INPUT_COORDINATES_FILE) %>% # Use read_tsv for tab-separated files
  dplyr::select(Probe_ID, CpG_chrm, CpG_beg, CpG_end) %>%
  filter(!is.na(CpG_chrm) & !is.na(CpG_beg) & !is.na(CpG_end)) %>%
  mutate(CpG_chrm = if_else(str_detect(CpG_chrm, "chr"), 
                           as.character(CpG_chrm), 
                           paste0("chr", CpG_chrm)))

# Create a GRanges object
coord_ranges <- GRanges(
  seqnames = Rle(input_data$CpG_chrm),
  ranges = IRanges(start = input_data$CpG_beg, end = input_data$CpG_end),
  name = input_data$Probe_ID
)

# --- Step 2: Load Gene Annotations and Find Overlaps ---
genes_db <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
overlaps <- findOverlaps(coord_ranges, genes_db)

# --- Step 3: Create and Save the Final Mapping File ---
probe_to_entrez_map <- data.frame(
  Probe_ID = coord_ranges$name[queryHits(overlaps)],
  EntrezID = genes_db$gene_id[subjectHits(overlaps)],
  stringsAsFactors = FALSE
) %>%
  distinct()

write_tsv(
  probe_to_entrez_map,
  file = OUTPUT_FILE
)

message(paste("Successfully wrote lookup table to:", OUTPUT_FILE))
