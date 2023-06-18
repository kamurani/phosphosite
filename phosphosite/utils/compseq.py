"""Util functions for sequence composition analysis."""
# %%
# Phosphosite
# Author: Cam İmran <c.mcmenamie@unsw.edu.au>
# License: MIT
# Project Website: https://github.com/kamurani/phosphosite
# Code Repository: https://github.com/kamurani/phosphosite

import os
import pandas as pd

from typing import Union
from pathlib import Path
from typing import Dict, List, Tuple, Union

from phosphosite import DATA_DIR 
DEFAULT_COMPSEQ_PATH = DATA_DIR / "compseq" / "scgr4_human.composition" 
AA_BACKGROUND_PATH = DATA_DIR / "compseq" / "aa_background.composition" 

DEFAULT_COMPSEQ_PATH


COMPSEQ_OUTPUT_SKIP = 26


def get_background_df(
    file_path: Union[str, Path] = AA_BACKGROUND_PATH,
) -> pd.DataFrame:
    """"Get background amino acid frequencies as a DataFrame.
    
    Parameters
    ----------
    file_path : Union[str, Path], optional
        The path to the background amino acid frequencies file, by default AA_BACKGROUND_PATH

    """
    df = parse_compseq_output(file_path)
    return df



"""Parse and load data from a `compseq` output file."""
def parse_compseq_output(
    file_path: Union[str, Path] = DEFAULT_COMPSEQ_PATH,

) -> pd.DataFrame: 
    """Parse compseq file.
    
    Output is a pandas DataFrame with the following rows:
    - `count` : The number of times the k-mer was observed in the input sequence(s).
    - `freq` : The frequency of the k-mer. 
    - `exp` : The expected frequency of the k-mer.
    - `freq/exp` : The ratio of the observed frequency to the expected frequency.

    Each column is a k-mer. 


    Parameters
    ----------
    file_path : Union[str, Path]
        The path to the `compseq` output file.

    Returns
    -------
    pd.DataFrame
        A pandas DataFrame containing the parsed data.
    """
    skip = COMPSEQ_OUTPUT_SKIP

    # Read in file
    df = pd.read_csv(
        file_path, 
        skiprows=skip,
        sep="\s+",
        header=None,
    )
    df.columns = ["word", "count", "freq", "exp", "freq/exp"]

    # Set 'word' to be index. 
    df = df.set_index("word")
    df.columns.name = "observations"
    df = df.T
    return df


def get_frequency_pair_dict(
    df: pd.DataFrame,
    amino_acids: str = "ACDEFGHIKLMNPQRSTVWY",
    metric: str = "freq",
) -> Dict[str, float]:
    """Get frequency pair dictionary.
    
    This assumes that word lengths are 3.
    We take the middle residue to be the generic residue ("X").
    For each pair (R1, R2), we sum the frequency of all R1-X-R2 triplets.

    Parameters
    ----------
    df : pd.DataFrame
        The DataFrame containing the compseq data.
    amino_acids : str, optional
        The amino acids to consider, by default "ACDEFGHIKLMNPQRSTVWY".
        If ``None`` is passed, all amino acids in the dataframe (possibly
        including non-standard amino acids) will be considered.
    metric : str, optional
        The metric to use, by default "freq".
        One of 
        - `count` : The number of times the k-mer was observed in the input sequence(s).
        - `freq` : The frequency of the k-mer.
        - `exp` : The expected frequency of the k-mer.
        - `freq/exp` : The ratio of the observed frequency to the expected frequency.

    
    Returns
    -------
    Dict[str, float]
        A dictionary containing the frequency pair data. 
        The keys are triplets (e.g. AXA), and the values are the sums of the specified metric.
         
        Note: the values are obtained from summing the individual triplet values. 


    """
    if amino_acids is None or len(amino_acids) == 0:
        amino_acids = get_all_amino_acids(df)

    freq_pair_dict = {}
    for aa_prev in amino_acids:
        for aa_next in amino_acids:
            # Each combination of pairs 

            # Get probability of AXA 
            triplets = [
                aa_prev + aa + aa_next
                for aa in amino_acids
            ]
            idx = aa_prev + "X" + aa_next
            freq_pair_dict[idx] = df[triplets].sum(axis=1)[metric] #.freq  # sum frequency for AXA 
    
    return freq_pair_dict

"""Get triplet dictionary for triplets that orient residues of interest in the middle."""
def get_triplet_dict(
    df: pd.DataFrame,
    amino_acids: Union[str, List] = "ACDEFGHIKLMNPQRSTVWY",
    metric: str = "freq",
    centre_residue: str = "STY",
) -> Dict[str, float]:
    """Get triplet dictionary.
    
    Parameters
    ----------
    df : pd.DataFrame
        The DataFrame containing the compseq data.
    amino_acids : str, optional
        The amino acids to consider, by default "ACDEFGHIKLMNPQRSTVWY".
        If ``None`` is passed, all amino acids in the dataframe (possibly
        including non-standard amino acids) will be considered.
    metric : str, optional
        The metric to use, by default "freq".
        One of 
        - `count` : The number of times the k-mer was observed in the input sequence(s).
        - `freq` : The frequency of the k-mer.
        - `exp` : The expected frequency of the k-mer.
        - `freq/exp` : The ratio of the observed frequency to the expected frequency.
    centre_residue : str, optional
        The residues to consider as the centre residue, by default "STY".
        If ``None`` is passed, all residues will be considered as the centre residue.
    
    Returns
    -------
    Dict[str, float]
        A dictionary containing the triplet data.          
        Note: the values are obtained from summing the individual triplet values. 
    """
    if amino_acids is None or len(amino_acids) == 0:
        amino_acids = get_all_amino_acids(df)

    triplets = []
    for res in centre_residue:
        for aa_prev in amino_acids:
            for aa_next in amino_acids:
                # Each combination of pairs 
                triplets.append(aa_prev + res + aa_next)
    
    stat_dict = df.to_dict(orient="index")
    metric_dict = stat_dict[metric]

    return {
        k: v 
        for k, v in metric_dict.items() 
        if k in triplets
    }

def get_all_amino_acids(
    df: pd.DataFrame,
    remove_other: bool = True,
) -> str: 
    """Get all unique residues in a triplet dataframe.
    
    Parameters
    ----------
    df : pd.DataFrame
        The DataFrame containing the compseq data.
    remove_other : bool, optional
        Whether to remove the "Other" column, by default True.
    
    Returns
    -------
    str
        A string containing all unique residues in the dataframe.
    """
    triplets = list(df.columns)
    if remove_other:
        triplets.remove("Other")
    amino_acids = []
    for i in [0, 1, 2]:
        amino_acids += list(set([t[i] for t in triplets]))
    amino_acids = sorted(list(set(amino_acids)))
    return "".join(amino_acids)