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
LOOK_UP_TABLE = APP_DIR / "lookup.tsv"


def annotate(lookup_df: pl.DataFrame, rdata_file_path: str, background_file_path: str, output_file: str):
    '''
    Annotes a background file with trait data from an RDS file.
    The result is a DF with a bitset for each probe ID, indicating which traits are associated with each EntrezID.
    '''
    # Load Data 
    probe_df = pl.read_csv(background_file_path, separator='\t', null_values='NA').drop_nulls(ENTEREZID_COL)
    raw_traits = read_rds(rdata_file_path)

    # --- Pre-compute EntrezID -> Bitmask Map in Python ---
    # First map entrez IDs to traits
    entrez_to_traits_map = defaultdict(list)
    for trait, eids in raw_traits.items():
        if eids:
            for eid_str in eids:
                try: entrez_to_traits_map[int(eid_str)].append(trait)
                except (ValueError, TypeError): continue

    trait_to_index_map = dict(zip(lookup_df["trait"], lookup_df["index"]))
    
    BITS_PER_CHUNK = 64
    num_chunks = ceil(lookup_df['index'].max() / BITS_PER_CHUNK)

    # Pre-calculate bitmasks for each EntrezID 
    entrez_bitmask_data = []
    for entrez_id, traits in entrez_to_traits_map.items():
        bitmask_array = [0] * num_chunks
        for trait in traits:
            if (bit_pos := trait_to_index_map.get(trait)) is not None:
                bitmask_array[bit_pos // BITS_PER_CHUNK] |= (1 << (bit_pos % BITS_PER_CHUNK))
        entrez_bitmask_data.append([entrez_id] + bitmask_array)

    # ---  Build Final DataFrame with Polars ---
    bitset_cols = [f"bitset_{i}" for i in range(num_chunks)]
    # Define the correct schema from the start
    schema = {ENTEREZID_COL: pl.Int64}
    schema.update({col: pl.UInt64 for col in bitset_cols})

    entrez_bitmask_df = pl.DataFrame(entrez_bitmask_data, schema=schema, orient="row")
   
    # Join, aggregate bitmasks with a bitwise OR, and create the final result in one chain.
    final_df = (
        probe_df.join(entrez_bitmask_df, on=ENTEREZID_COL, how="inner")
        .group_by(PROBE_ID_COL)
        .agg([pl.col(col).bitwise_or() for col in bitset_cols])
    ).sort(PROBE_ID_COL)

    final_df.write_csv(output_file, separator='\t', null_value='NA')


def create_lookup_table(lookup_table:str, rdata_file:str):
    if Path(lookup_table).exists():
        return pl.read_csv(lookup_table, separator='\t', null_values='NA')

    content = read_rds(rdata_file)
    keys = sorted(content.keys())
    values = [*range(len(keys))]
    df = pl.DataFrame({'trait': keys, 'index': values})
    df.write_csv(lookup_table, separator='\t', null_value='NA')
    return df


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
    
    lookup_df = create_lookup_table(LOOK_UP_TABLE, args.rdata_file)
    annotate(
        lookup_df=lookup_df,
        rdata_file_path=args.rdata_file,
        background_file_path=args.input_file,
        output_file=args.output_file
    )


if __name__ == "__main__":
    main()