"""Process bubbles for a proteome."""
# TODO:
# - add click command line interface for supplying a text file of uniprot ids


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
from phosphosite.bubble import process
        



phosphosite_df = None

# Read in uniprot_ids 
list_path = Path("./uniprot_ids.txt")
with open(list_path, "r") as f:
    uniprot_ids = f.read().splitlines()

"""Save path"""
result_path = Path("./results/bubbles.h5")


adjacency_range = 2
radius = 8.0 #6.0
kwargs = {
    "protein_ids": uniprot_ids,
    "phosphosite_df": phosphosite_df,
    "adjacency_range": adjacency_range,
    "radius": radius,
    "verbose": False,
    "filepath": result_path,
    "overwrite": True,
}
result = process(**kwargs)