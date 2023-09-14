"""Util functions for reading and writing data."""

from pathlib import Path
from typing import Dict, List, Union, Any

from Bio import SeqIO


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