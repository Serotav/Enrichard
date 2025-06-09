import pandas as pd
import numpy as np
import os
import rpy2.robjects as robjects
from collections import defaultdict
from sys import argv, stderr

df = pd.read_csv("MSA_annotation/MSA.hg38.manifest.gencode.v41.tsv",sep ='\t')

if 'info' in argv: print(f'{df.shape=}',file=stderr)
df = df[df["genesUniq"].notna()]

intermediate = open('background/symbols_to_entrezid.txt').read().strip().split('\n')
mapping = {symbol: entreid for line in intermediate for (symbol, entreid) in [line.split(" ")]}

enterzid = []
for gene in df["genesUniq"]:
    enterzid.append(';'.join([mapping[symbol] for symbol in gene.split(';') if symbol in mapping]))

df["ENTREZID"] = enterzid
if 'info' in argv: print(f'After enterzid {df.shape=}',file=stderr)

def load_rdata(data:str):
    if 'info' in argv: print(f'Partsing: {data}',file=stderr)
    robjects.r['load'](data)
    last_obj_name = str(robjects.r('ls()')[-1])
    target = robjects.r[last_obj_name]
    robjects.r(f'rm({last_obj_name})')
    return {key: list(target.rx2(key)) for key in target.names}

data = []
for file_name in sorted(os.popen('find rdata -type f').read().split()):
    data.append((file_name,load_rdata(file_name)))

for file in data:
    file_name, content = file[0], file[1]
    eid_to_trait = defaultdict(lambda: list())
    for trait in content.keys():
        for eid in content[trait]:
            eid_to_trait[eid].append(trait)
    v = [np.nan] * df.shape[0]
    for i,eid in enumerate(df["ENTREZID"]):
        if eid in eid_to_trait:
            v[i] =';'.join(eid_to_trait[eid])
    df[file_name.split('/')[-1]] = v
    if 'info' in argv: print(f'{file_name.split('/')[-1]} {df[file_name.split('/')[-1]].count()}',file=stderr)

for col in df.columns:
    if pd.api.types.is_numeric_dtype(df[col]) and not df[col].isna().any():
        df[col] = df[col].astype(int)  


df.to_csv("background/MSA_parsed.csv",index=False, sep='\t')
