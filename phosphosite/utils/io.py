"""Util functions for writing and reading files."""
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

COMPSEQ_OUTPUT_SKIP = 26


"""Parse and load data from a `compseq` output file."""
def parse_compseq(
    file_path: Union[str, Path], 

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