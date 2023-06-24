"""Phosphosite package."""

import os
from pathlib import Path

PROJECT_ROOT_DIR = Path(__file__).resolve().parent.parent
PHOSPHOSITE_DIR = PROJECT_ROOT_DIR / 'phosphosite'
DATA_DIR        = PROJECT_ROOT_DIR / 'data'
DATASET_DIR     = PROJECT_ROOT_DIR / 'datasets'
PSP_DATASET_DIR = DATASET_DIR / 'psp'

STRUCTURE_DIR   = DATA_DIR / 'structures'
CIF_DIR         = STRUCTURE_DIR / 'alphafold' / 'cif'
PAE_DIR         = STRUCTURE_DIR / 'alphafold' / 'pae'
PDB_DIR         = STRUCTURE_DIR / 'alphafold' / 'pdb' # PDB format structures from AF database.
SEQUENCE_DIR    = DATA_DIR / 'sequence'

NOTEBOOK_DIR    = PROJECT_ROOT_DIR / 'notebooks'

CONFIG_FILE_PATH = NOTEBOOK_DIR / "config.yml"

import yaml 
with open(CONFIG_FILE_PATH, 'r') as ymlfile:
    cfg = yaml.load(ymlfile, Loader=yaml.FullLoader)

# Get the path to the DATASET_DIR
USER_DATASET_DIR = Path(cfg['DATASET_DIR']).expanduser()
USER_STRUCTURE_DIR = Path(cfg['STRUCTURE_DIR']).expanduser()

"""PhosphoSitePlus datasets""" 
PSP_DIR = USER_DATASET_DIR / 'PSP' 
assert PSP_DIR.is_dir()

P_PATH          = PSP_DIR / "Phosphorylation_site_dataset"  # Phosphorylation sites dataset
REG_PATH        = PSP_DIR / "Regulatory_sites"              # Regulatory sites dataset
PTM_SEQ_PATH    = PSP_DIR / "Phosphosite_PTM_seq.fasta"     # FASTA formatted sequences of all substrate proteins.  Lowercase letters indicate the phosphorylation sites.

# Assert that the paths are existing files
assert P_PATH.is_file()
assert REG_PATH.is_file()
assert PTM_SEQ_PATH.is_file()

"""AlphaFold downloaded structures"""
AF_HUMAN_CIF = USER_STRUCTURE_DIR


AF_VERSION = 3
DEFAULT_RADIUS = 6.0

if __name__ == '__main__':
    print(f"PROJECT_ROOT_DIR: {PROJECT_ROOT_DIR}")
    print(f"DATA_DIR: {DATA_DIR}")