"""Phosphosite package."""

import os
from pathlib import Path

PROJECT_ROOT_DIR = Path(__file__).resolve().parent.parent
PHOSPHOSITE_DIR = PROJECT_ROOT_DIR / 'phosphosite'
DATA_DIR        = PROJECT_ROOT_DIR / 'data'

STRUCTURE_DIR   = DATA_DIR / 'structures'
CIF_DIR         = STRUCTURE_DIR / 'alphafold' / 'cif'
PAE_DIR         = STRUCTURE_DIR / 'alphafold' / 'pae'
PDB_DIR         = STRUCTURE_DIR / 'alphafold' / 'pdb' # PDB format structures from AF database.
SEQUENCE_DIR    = DATA_DIR / 'sequence'


if __name__ == '__main__':
    print(f"PROJECT_ROOT_DIR: {PROJECT_ROOT_DIR}")
    print(f"DATA_DIR: {DATA_DIR}")