import polars as pl
from pathlib import Path
from collections import defaultdict
from sys import stderr
from rds2py import read_rds
import argparse
from math import ceil
import pathlib
"""
- Iterates through each RDS file in `rdata/`, reading via rds2py.
  - Each RDS is a dict: trait → [Entrez IDs].
  - Converts it to Entrez ID → traits mapping, then annotates the DataFrame.
"""
ENTEREZID_COL = "EntrezID"
PROBE_ID_COL = "Probe_ID"
APP_DIR = pathlib.Path(__file__).parent


def annotate(rdata_file_path: str, background_file_path: str, output_file: str):
    '''
    Annotes a background file with trait data from an RDS file.
    The result is a DF with a a col for each trait.
    '''
    # Load Data 
    probe_to_entrez_df = pl.read_csv(background_file_path, separator='\t', null_values='NA').drop_nulls(ENTEREZID_COL)
    raw_traits = read_rds(rdata_file_path)

    # Mapping Entrez IDs to traits
    trait_to_entrez = {}
    for trait, entrez_ids in raw_traits.items():
        trait_to_entrez[trait] = [int(eid) for eid in entrez_ids]

    trait_df_intermediate = pl.DataFrame([trait_to_entrez])


    trait_df = trait_df_intermediate.unpivot(
    index=[], on=list(trait_to_entrez.keys())
    ).rename(
        {'variable': 'Trait', 'value': 'EntrezID'}
    ).explode(
        'EntrezID'
    )

    
    merged_df = probe_to_entrez_df.join(trait_df, on='EntrezID', how='inner').with_columns(
        pl.lit(1).cast(pl.UInt8).alias("has_trait")
    )
    
    final_df = merged_df.pivot(
        index="Probe_ID",
        on="Trait",
        values="has_trait",
        aggregate_function="max"
    ).fill_null(0)

    final_df.write_parquet(output_file, compression='lz4')


def parse_args():
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
        'rdata_file', 
        help='Directory containing RDS files'
    )


    
    return parser.parse_args()


def main():
    args = parse_args()

    if not pathlib.Path(args.input_file).exists():
        print(f"Input file {args.input_file} does not exist.", file=stderr)
        return
    if not pathlib.Path(args.rdata_file).exists():
        print(f"RDS directory {args.rdata_dir} does not exist or is not a directory.", file=stderr)
        return
    
    annotate(
        rdata_file_path=args.rdata_file,
        background_file_path=args.input_file,
        output_file=args.output_file
    )


if __name__ == "__main__":
    main()