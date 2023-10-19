import pandas as pd 
import numpy as np

from pathlib import Path
from phosphosite import UNIPROT_SEQUENCE_PATH
from phosphosite.utils import load_seq_dict_from_file


"""Sequence dictionary"""
sequence_dict = load_seq_dict_from_file(UNIPROT_SEQUENCE_PATH, key_format="uniprot_id")


