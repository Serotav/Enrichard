import pandas as pd
import numpy as np
import pathlib
from collections import defaultdict
from sys import argv, stderr
from rds2py import read_rds
import argparse
"""
- Iterates through each RDS file in `rdata/`, reading via rds2py.
  - Each RDS is a dict: trait → [Entrez IDs].
  - Converts it to Entrez ID → traits mapping, then annotates the DataFrame.
"""
ENTEREZID_COL = "EntrezID"

def annotate(input_file:str, output_file:str, rdata_dir:str):
    # Loading the df we want to annotate
    if input_file.split('.')[-1] == 'csv':
        background_df = pd.read_csv(input_file)
    elif input_file.split('.')[-1] == 'tsv':
        background_df = pd.read_csv(input_file, sep='\t')
    else:
        print(f"Unsupported file format: {input_file}", file=stderr)
        return

    # Cast EntrezID column to numeric this should not be necessary, but R is not trustworthy
    background_df[ENTEREZID_COL] = pd.to_numeric(background_df[ENTEREZID_COL], errors='coerce')
    
    print(f'Parsing {input_file} {background_df.columns}',file=stderr)

    # We can load the data from bioinf.wehi.edu.au
    # And parse the data, for each file (trait) we will create a new column in the df
    for file_name in sorted(str(f) for f in pathlib.Path(rdata_dir).rglob("*") if f.is_file()):
        # Eeach file maps trait -> entrezid, but we need to map entrezid -> trait(s)
        # So we will create a dict where the key is the entrezid and the value is a list of traits
        print(f'Processing {file_name}', file=stderr)
        content = read_rds(file_name)
        eid_to_traits = defaultdict(lambda: list())
        
        for trait in content.keys():
            for eid in content[trait]:
                eid_to_traits[int(eid)].append(trait)
        
        def join_traits(trait_list):
            return ';'.join(trait_list) if trait_list else np.nan

        # Map the EntrezIDs to their traits using the dictionary
        col_name = file_name.split('/')[-1]
        background_df[col_name] = background_df[ENTEREZID_COL].map(eid_to_traits).apply(join_traits)

        if 'info' in argv: print(f'{file_name.split('/')[-1]} {background_df[file_name.split('/')[-1]].count()}',file=stderr)

    # Casting bullshit
    for col in background_df.columns:
        if pd.api.types.is_numeric_dtype(background_df[col]) and not background_df[col].isna().any():
            background_df[col] = background_df[col].astype(int)  

    print(f'Final shape {background_df.shape} location {output_file}',file=stderr)
    background_df.to_csv(output_file,index=False, sep='\t')


def main():
    parser = argparse.ArgumentParser(
        description="Annotate a file with trait data from RDS files.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    parser.add_argument(
        'input_file',
        help='Input file to annotate (CSV or TSV format)'
    )
    
    parser.add_argument(
        'output_file', 
        help='Output file name'
    )

    parser.add_argument(
        'rdata_dir', 
        help='Output file name'
    )

    args = parser.parse_args()

    annotate(args.input_file, args.output_file, args.rdata_dir)


if __name__ == "__main__":
    main()