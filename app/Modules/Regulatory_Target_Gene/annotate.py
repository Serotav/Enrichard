import polars as pl
from pathlib import Path
from sys import stderr
from rds2py import read_rds
import argparse
import pathlib
import time
"""
This script reads an RDS file containing trait data and a background file
with probe IDs and Entrez IDs. It then annotates the background file with
the trait data, creating a DataFrame where each row corresponds to a probe ID
and each column corresponds to a trait. The output is saved as a Parquet file.
"""
ENTEREZID_COL = "EntrezID"
PROBE_ID_COL = "Probe_ID"
APP_DIR = pathlib.Path(__file__).parent


def annotate(rdata_file_path: Path, background_file_path: Path, output_file: Path) -> pl.DataFrame:
    '''
    Annotes a background file with trait data from an RDS file.
    The result is a DF with a Probe_Id followed by col for each trait.
    
    NB: we need to handle cases where the same cg probes maps to multiple Entrez IDs.
        and where the same enterez ID maps to multiple cg probes.
    '''
    # Load Data 
    # Probe_ID | EntrezID
    # id1      | 1234
    # id2      | 5678
    probe_to_entrez_df = pl.read_csv(background_file_path, separator='\t', null_values='NA').drop_nulls(ENTEREZID_COL)
    raw_traits = read_rds(str(rdata_file_path))

    # Mapping traits -> [Entrez IDs] Enterzids are integers
    trait_to_entrez = {}
    for trait, entrez_ids in raw_traits.items():
        trait_to_entrez[trait] = [int(eid) for eid in entrez_ids]

    # Convert to DataFrame
    # Trait1 | Trait2 
    # [id1, id2]    | [id4, id5, id1]
    trait_df_intermediate = pl.DataFrame([trait_to_entrez])

    # Long format:
    # Trait      | EntrezID
    # Trait1     | id1
    # Trait1     | id2
    # Trait2     | id4
    trait_df = trait_df_intermediate.unpivot(
        index=[], on=list(trait_to_entrez.keys())
    ).rename(
        {'variable': 'Trait', 'value': 'EntrezID'}
    ).explode('EntrezID')
    
    # adding column Probe_ID
    # Trait      | EntrezID | Probe_ID
    # Trait1     | id1      | cg1   
    # Trait1     | id2      | cg2   
    # Trait2     | id4      | cg3   
    merged_df = probe_to_entrez_df.join(trait_df, on='EntrezID', how='inner').with_columns(
        pl.lit(1).cast(pl.UInt8).alias("has_trait")
    )
    
    # Pivoting the traits to columns
    # Probe_ID  | Trait1 | Trait2 | Trait3
    # cg1       | 1      | 0      | 0
    # cg2       | 1      | 1      | 0
    final_df = merged_df.pivot(
        index="Probe_ID",
        on="Trait",
        values="has_trait",
        aggregate_function="max"
    ).fill_null(0)

    final_df.write_parquet(output_file, compression='lz4')
    return final_df

def cache_annotated_df(annotated_df: pl.DataFrame, cached_file: Path) -> None:
    # Unpivoting to one trait per row, summing how many probes are associated with each trait.
    # Trait      | totals
    # Trait1     | 69
    # Trait2     | 420
    contingency_table_df = annotated_df.select(
        [pl.col(col).sum().alias(col)for col in annotated_df.columns if col != PROBE_ID_COL]
    ).unpivot(index=[], variable_name="Trait", value_name="totals")

    contingency_table_df.write_parquet(cached_file, compression='lz4')


def parse_args()-> argparse.Namespace:
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
        help='The RDS file containing traits -> EnterezIds'
    )
    parser.add_argument(
        'cached_file', 
        help='Cached output file'
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
    
    start = time.time()
    anotated_df = annotate(
        rdata_file_path= Path(args.rdata_file),
        background_file_path=Path(args.input_file),
        output_file=Path(args.output_file)
    )
    print(f"Annotation completed in {time.time() - start:} seconds.")
    
    start = time.time()
    cache_annotated_df(anotated_df, Path(args.cached_file))
    print(f"Caching completed in {time.time() - start:} seconds.")

if __name__ == "__main__":
    main()