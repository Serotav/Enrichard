import pandas as pd
import numpy as np
import pathlib
from collections import defaultdict
from sys import argv, stderr
from rds2py import read_rds

"""
TL;DR: we take genesUniq and map them to traits.
- Loads a Gencode manifest TSV, filters rows with non-null gene symbols.
- Builds a gene-symbol → Entrez ID mapping from a lookup file.
- Maps each rows semicolon-separated gene symbols to Entrez IDs.
- Iterates through each RDS file in `rdata/`, reading via rds2py.
  - Each RDS is a dict: trait → [Entrez IDs].
  - Converts it to Entrez ID → traits mapping, then annotates the DataFrame.
- Converts fully numeric columns to integers.
- Outputs the final annotated table to `background/MSA_parsed.csv`.
"""

if __name__ == "__main__":
    # Loading the df we want to annotate
    background_df = pd.read_csv("MSA_annotation/MSA.hg38.manifest.gencode.v41.tsv",sep ='\t')
    background_df = background_df[background_df["genesUniq"].notna()]
    if 'info' in argv: print(f'{background_df.shape=}',file=stderr)

    # Create the mapping dict, here the mapping is symbol/ENSEMBL ->  entrezid
    mapping_file = open('background/symbols_to_entrezid.txt').read().strip().split('\n')
    mapping = {symbol: entreid for line in mapping_file for (symbol, entreid) in [line.split(" ")]}

    # Each cpg site might have might have multiple "genesUniq" entries, so we need to map each of them to entrezid
    # Some might no have a mapping, TODO: drop those
    enterzid = []
    for gene in background_df["genesUniq"]:
        enterzid.append(';'.join([mapping[symbol] for symbol in gene.split(';') if symbol in mapping]))

    # Once mapped, we can add the entrezid to the df. The mapping is complete
    background_df["ENTREZID"] = enterzid
    if 'info' in argv: print(f'After enterzid {background_df.shape=}',file=stderr)

    # Now we can load the data from bioinf.wehi.edu.au
    # And parse the data, for each file (trait) we will create a new column in the df
    for file_name in sorted(str(f) for f in pathlib.Path("rdata").rglob("*") if f.is_file()):
        # Eeach file maps trait -> entrezid, but we need to map entrezid -> trait(s)
        # So we will create a dict where the key is the entrezid and the value is a list of traits
        print(f'Processing {file_name}', file=stderr)
        content = read_rds(file_name)
        eid_to_traits = defaultdict(lambda: list())
        
        for trait in content.keys():
            for eid in content[trait]:
                eid_to_traits[eid].append(trait)
        
        # Now we can create a new column in the background_df for each file
        # The column will contain the traits for the entrezid in that row
        new_col = [np.nan] * background_df.shape[0]
        for i,eid in enumerate(background_df["ENTREZID"]):
            if eid in eid_to_traits:
                new_col[i] =';'.join(eid_to_traits[eid])
        
        background_df[file_name.split('/')[-1]] = new_col

        if 'info' in argv: print(f'{file_name.split('/')[-1]} {background_df[file_name.split('/')[-1]].count()}',file=stderr)

    # Casting bullshit
    for col in background_df.columns:
        if pd.api.types.is_numeric_dtype(background_df[col]) and not background_df[col].isna().any():
            background_df[col] = background_df[col].astype(int)  

    if 'info' in argv: print(f'Final shape {background_df.shape}', file=stderr)
    background_df.to_csv("background/MSA_parsed.csv",index=False, sep='\t')
