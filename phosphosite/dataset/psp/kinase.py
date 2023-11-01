"""Kinase substrate dataset from PSP."""
# %%
# Phosphosite
# Author: Cam İmran <c.mcmenamie@unsw.edu.au>
# License: MIT
# Project Website: https://github.com/kamurani/phosphosite
# Code Repository: https://github.com/kamurani/phosphosite
import pandas as pd 
import numpy as np
import re

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from typing import Union, List, Tuple, Dict
from pathlib import Path

from phosphosite import PSP_DATASET_DIR
from pathlib import Path


# Kinase substrate dataset
filepath = PSP_DATASET_DIR / "Kinase_Substrate_Dataset.gz"
kinase_substrate_dataset = pd.read_csv(filepath, sep="\t", compression="gzip")

