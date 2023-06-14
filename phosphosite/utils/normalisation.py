"""Normalise data."""
# %%
# Phosphosite
# Author: Cam İmran <c.mcmenamie@unsw.edu.au>
# License: MIT
# Project Website: https://github.com/kamurani/phosphosite
# Code Repository: https://github.com/kamurani/phosphosite

from phosphosite.utils.compseq import parse_compseq_output

import pandas as pd 

freq_df = parse_compseq_output()

"""Get background frequency of amino acids as matrix."""
def get_frequency_matrix(
    df: pd.DataFrame = freq_df,
    centre_residue: str = "S",
    metric: str = "freq", # "freq/exp"

) -> pd.DataFrame: 
    """Generates a matrix representing triplet background frequency values.
    
    Each cell represents a triplet frequency from the triplet ROW-centre_residue-COLUMN.


    """
    # Get all triplets
    triplets = df.columns

    # Get all residues
    residues = set()
    for triplet in triplets:
        residues.add(triplet[0])
        residues.add(triplet[2])
    residues = sorted(list(residues))

    # Generate matrix
    matrix = pd.DataFrame(index=residues, columns=residues)
    for triplet in triplets:
        matrix.loc[triplet[0], triplet[2]] = df.loc[metric, triplet]
    
    # Fill in missing values
    for i in residues:
        for j in residues:
            if i == centre_residue or j == centre_residue:
                continue
            if matrix.loc[i, j] is None:
                matrix.loc[i, j] = matrix.loc[j, i]
    return matrix