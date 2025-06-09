library(org.Hs.eg.db)
library(AnnotationDbi)
library(tidyverse)


all_ensembl <- keys(org.Hs.eg.db, keytype = "ENSEMBL")
all_symbols <- keys(org.Hs.eg.db, keytype = "SYMBOL")


dioporco <-mapIds(org.Hs.eg.db, keys = all_ensembl, column = "ENTREZID", keytype = "ENSEMBL",multiVals = "first")
diocane <-mapIds(org.Hs.eg.db, keys = all_symbols, column = "ENTREZID", keytype = "SYMBOL",multiVals = "first")

porcamadonna <-c(dioporco,diocane) 
df <- data.frame(name = names(porcamadonna), value = porcamadonna, row.names = NULL)
write.table(df, file = "background/symbols_to_entrezid.txt", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)
