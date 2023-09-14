import pandas as pd 
import numpy as np 
import h5py

from tqdm import tqdm
from pathlib import Path
from typing import Dict, List, Union, Any


from phosphosite import PER_RESIDUE_EMBEDDING_PATH 


"""Retrieve a single embedding from h5 file."""
def get_embedding(
    uniprot_id: str, 
    embeddings_path: Path = PER_RESIDUE_EMBEDDING_PATH,
    dtype: str = "float32",

    average: bool = False,
) -> Any:
    with h5py.File(embeddings_path, "r") as f: 
        try:
            embedding = f[uniprot_id][()]
        except KeyError:
            raise ValueError(f"{uniprot_id}: No embedding data.")

        if average: embedding = np.mean(embedding, axis=0)

        if dtype is not None:
            embedding = embedding.astype(dtype)
        return embedding
            



