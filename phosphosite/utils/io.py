"""Util functions for reading and writing data."""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Union, Any

from Bio import SeqIO
# Download fasta files if they don't exist.

import requests
import Bio.SeqIO as SeqIO
import time 

from pathlib import Path
from typing import List, Dict, Union, Tuple
from phosphosite import UNIPROT_DIR
from phosphosite.utils.amino import standard_amino_acids, amino_acid_names




def download_fasta(uniprot_ids: List[str], sequence_dir: Path, overwrite: bool = False):
    """Download fasta files from uniprot."""
    for uniprot_id in uniprot_ids:
        filename = "{uniprot_id}.fasta" 
        fp = sequence_dir / filename.format(uniprot_id=uniprot_id)

        if not fp.exists() or overwrite:
            url = f"https://www.uniprot.org/uniprot/{uniprot_id}.fasta"
            r = requests.get(url)
            r.raise_for_status()

            with open(fp, "w") as f:
                f.write(r.text)

# Load

def read_sequence_from_file(uniprot_id: str, sequence_dir: Path = UNIPROT_DIR):
    """Returns sequence from fasta file."""
    fasta_file = sequence_dir / f"{uniprot_id}.fasta"
    assert fasta_file.exists(), f"{fasta_file} does not exist."

    with open(fasta_file, "r") as f:
        record = SeqIO.read(f, "fasta")
        return str(record.seq)

def load_seq_dict_from_dir(
    sequence_dir: Path,
) -> Dict[str, str]:
    """Create dictionary for sequences from fasta files.
    
    Keys are uniprot ids, values are sequences.

    Parameters
    ----------
    sequence_dir : Path, optional
        Directory containing fasta files, by default None
    
    """
    if not sequence_dir.is_dir():
        raise ValueError(f"{sequence_dir} is not a directory.")
    # Load in each fasta file in the directory
    seq_dict = {}
    for fasta_file in UNIPROT_DIR.glob("*.fasta"):
        with open(fasta_file, "r") as f:
            record = SeqIO.read(f, "fasta")
            seq_dict[fasta_file.stem] = str(record.seq)
        
    return seq_dict

def load_seq_dict_from_file(
    sequence_file: Path,
    key_format: str = "uniprot_id",
) -> Dict[str, str]:
    """Create dictionary for sequences from single fasta file.
    
    Keys are uniprot ids, values are sequences.

    Parameters
    ----------
    sequence_file : Path, optional
        Fasta file, by default None
    
    """
    if not sequence_file.is_file():
        raise ValueError(f"{sequence_file} is not a file.")
    # Load in each fasta file in the directory
    seq_dict = {}
    with open(sequence_file, "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            if key_format == "uniprot_id":
                key = record.id.split("|")[1]
            else: 
                key = record.id

            seq_dict[key] = str(record.seq)
        
    return seq_dict
    

"""Parse peptide sequence"""
def parse_peptide_seq(peptide_seq: str) -> str:
    """Parse peptide sequence and return list of modifications and uniprot_id."""
    split = peptide_seq.split("_")
    uniprot_id = split[0]
    sites = split[1:]

    return uniprot_id, sites


"""Generate hash id using name and current time."""
def generate_hash_id(
    name: str, 
    upper: bool = True, 
    first_n_chars: int = 6,
    negative: bool = False,
):

    # current date and time including milliseconds
    now = time.time_ns()
    id = hash(name + str(now))

    if not negative: id = abs(id)
    # hex value of the id
    id = hex(id)
    if upper: id = id.upper()
    return id[0:first_n_chars]


import torch
import json

def save_index_dict(
    index_dict: Dict[str, Dict[str, torch.Tensor]],
    filepath: Path,
    overwrite: bool = False,
) -> None:
    """ 
    Save a dictionary of indexes to a json file. 

    NOTE: convert Tensors to lists before saving. 
    """

    # Create a copy of the index_dict with lists instead of tensors

    idx_dict = index_dict.copy()

    for key, value in idx_dict.items():
        idx_dict[key]["idx"] = value["idx"].tolist()
        idx_dict[key]["y"] = value["y"].tolist()


    if filepath.exists() and not overwrite:
        raise ValueError(f"File already exists at {filepath}")
    with open(filepath, "w") as f:
        json.dump(idx_dict, f)    

def load_index_dict(
    filepath: Path,
) -> Dict[str, Dict[str, torch.Tensor]]:
    """ 
    Load a dictionary of indexes from a json file. 

    NOTE: convert lists to Tensors after loading. 
    """
    with open(filepath, "r") as f:
        index_dict = json.load(f)    

    for key, value in index_dict.items():
        index_dict[key]["idx"] = torch.tensor(value["idx"])
        index_dict[key]["y"] = torch.tensor(value["y"])

    return index_dict