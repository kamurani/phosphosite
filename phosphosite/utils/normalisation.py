"""Normalise data."""
# %%
# Phosphosite
# Author: Cam İmran <c.mcmenamie@unsw.edu.au>
# License: MIT
# Project Website: https://github.com/kamurani/phosphosite
# Code Repository: https://github.com/kamurani/phosphosite

from phosphosite.utils.compseq import parse_compseq_output, get_frequency_pair_dict
from phosphosite.utils import standard_residues, aa1to3
import pandas as pd 


"""3D normalisation."""
# TODO



"""Get background frequency of amino acids as matrix."""
def get_frequency_matrix(
    df: pd.DataFrame,
    centre_residue: str = "S", # "X"
    metric: str = "freq", # "freq/exp"
    columns: str = standard_residues,
    letter_code: int = 3, # or 1

) -> pd.DataFrame: 
    """Generates a 2D matrix representing triplet background frequency values.
    
    Each cell represents a triplet frequency from the triplet ROW-centre_residue-COLUMN.

    Parameters
    ----------
    df : pd.DataFrame
        Triplet frequency dataframe. 
    centre_residue : str, optional
        The centre residue, by default "S".
    metric : str, optional
        The metric to use, by default "freq".
    columns : str, optional
        The residues to include, by default standard_residues. 
        Columns and rows will be ordered by this list.


    """
    # Get all triplets
    triplets = df.columns

    # Get all residues
    residues = set()
    for triplet in triplets:
        residues.add(triplet[0])
        residues.add(triplet[2])
    residues = sorted(list(residues))
    if columns is None: 
        columns = residues
    elif isinstance(columns, str):
        columns = list(columns)
    
    if centre_residue == "X":
        freq_dict = get_frequency_pair_dict(df, metric=metric) # for referencing like 'AXA'
    else: 
        freq_dict = df.loc[metric] # behaves like dict; for referencing like 'ASA'

    

    # Generate matrix
    matrix = pd.DataFrame(index=columns, columns=columns)
    for prev_res in columns: 
        for next_res in columns:
        
            #prev_res = triplet[0]   # ROW
            #next_res = triplet[2]   # COLUMN
            triplet_idx = prev_res + centre_residue + next_res 
            matrix.loc[prev_res, next_res] = freq_dict[triplet_idx]
    
    
    # Sort columns and rows

    if letter_code == 3: 
        # rename columns to 3 letter code
        matrix.columns = [aa1to3[aa] for aa in matrix.columns]
        matrix.index = [aa1to3[aa] for aa in matrix.index]
    
    # Sort columns and rows
    matrix = matrix.sort_index(axis=0)
    matrix = matrix.sort_index(axis=1)


    
    matrix.index.name = "prev"
    matrix.columns.name = "next"
    return matrix