# Here we are building a look up table symbol/ENSEMBL ->  entrezid
# Since the background could change, but we don't want to run this multiple times
# We will dump all the symbols and entrezids to a file regardless of them being in the background
# TODO: rewrite this in python
library(org.Hs.eg.db)
library(AnnotationDbi)
library(tidyverse)


# Get all ENSEMBL gene IDs
all_ensembl_ids <- keys(org.Hs.eg.db, keytype = "ENSEMBL")
#  Get all l Gene Symbols 
all_gene_symbols <- keys(org.Hs.eg.db, keytype = "SYMBOL")

# Look up the Entrez ID for each ENSEMBL ID. Returns ENSEMBL -> ENTREZID
ensembl_to_entrez_map <-mapIds(org.Hs.eg.db, keys = all_ensembl, column = "ENTREZID", keytype = "ENSEMBL",multiVals = "first")
# Look up the Entrez ID for each SYMBOL. Returns SYMBOL -> ENTREZID
symbol_to_entrez_map <-mapIds(org.Hs.eg.db, keys = all_symbols, column = "ENTREZID", keytype = "SYMBOL",multiVals = "first")

# Paste them together and write to a file
combined_id_map <- c(ensembl_to_entrez_map, symbol_to_entrez_map)
df <- data.frame(name = names(combined_id_map), value = combined_id_map, row.names = NULL)
write.table(df, file = "background/symbols_to_entrezid.txt", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)
