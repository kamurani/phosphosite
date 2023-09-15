"""Process bubbles for a proteome."""
# TODO:
# - add click command line interface for supplying a text file of uniprot ids


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
        





import click as ck 

@ck.command()
@ck.option("--uniprot-ids", "-u", help="Path to uniprot ids text file.", required=True,
           type=ck.Path(exists=True, dir_okay=False, path_type=Path)) 
@ck.option("--output", "-o", help="Path to output file.", required=True,
           type=ck.Path(dir_okay=True))
@ck.option("--overwrite/--no-overwrite", "-w/-nw", help="Overwrite existing file.", default=False)

def main(
    uniprot_ids, 
    output,
    overwrite,
):
    """Process bubbles for a proteome."""
    output = Path(output)
    # If output is directory, then save to default filename
    if output.is_dir():
        output = output / "all_bubbles.h5"
    
    # Check if output file exists and is not empty
    if output.exists() and output.stat().st_size > 0: 
        # Load in dataframe
        dataframe = pd.read_hdf(output)
        ck.echo(f"Loaded existing data from {output} with {len(dataframe['protein_id'].unique())} uniprot ids.")
        if overwrite:
            ck.echo("Overwriting existing data.")
        else:
            ck.echo("Appending to existing data.") 
    # Read in uniprot_ids 
    with open(uniprot_ids, "r") as f:
        uniprot_ids = f.read().splitlines()
    phosphosite_df = None

    adjacency_range = 2
    radius = 8.0 #6.0
    kwargs = {
        "protein_ids": uniprot_ids,
        "result_df": dataframe, 
        "overwrite": overwrite,
        "phosphosite_df": phosphosite_df,
        "adjacency_range": adjacency_range,
        "radius": radius,
        "verbose": False,
        "filepath": output,
        "overwrite": True,
    }
    result = process(**kwargs)


if __name__ == "__main__":
    main()