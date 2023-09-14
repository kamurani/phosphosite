"""Load saved bubble data."""


from pathlib import Path 
import pandas as pd 
import numpy as np
import os 
import re 
import gzip 
import shutil
import Bio.PDB.MMCIF2Dict
from typing import Union, List, Tuple, Dict, Optional
from pathlib import Path
from tqdm import tqdm

pd.options.mode.chained_assignment = None  # default='warn'

from phosphosite.utils import aa1to3, aa3to1
from phosphosite import NOTEBOOK_DIR 

verbose = False

result_path = NOTEBOOK_DIR / "bubble" / "results" / "all_bubbles.h5" 

## Load results
result_df = pd.read_hdf(result_path, key="data")
if verbose: print(f"Loaded {len(result_df.protein_id.unique())} unique proteins from {result_path}.")

## Annotate with sequence distance
result_df["seq_dist"] = result_df["nn_pos"] - result_df["pos"]

## Annotate with triplet 
result_df["triplet"] = result_df["prev"] + result_df["res"] + result_df["next"]

### Annotate known phosphosites 

from phosphosite.dataset.psp import phosphorylation
df = phosphorylation[["ACC_ID", "MOD_RSD"]].drop_duplicates()
df.rename(columns={"ACC_ID": "protein_id"}, inplace=True)
df["modification"] = df["MOD_RSD"].str.split("-").str[1]
df["res"] = df["MOD_RSD"].str.split("-").str[0].str[0]
df["pos"] = df["MOD_RSD"].str.split("-").str[0].str[1:].astype(int)
df = df[["protein_id", "res", "pos", "modification"]]

# Create boolean column for whether a residue is phosphorylated
phos_series = result_df.merge(df, on=["protein_id", "res", "pos"], how="left").modification.notna() #.astype(int) 
result_df["phosphosite"] = phos_series.values 

if verbose: print(f"Annotated {phos_series.sum()} phosphorylated residues.")
