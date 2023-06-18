"""Loading data for phosphosite analysis."""


from phosphosite import DATA_DIR

import pandas as pd
import numpy as np
from pathlib import Path

def load_saved_motif_df(
    to_process: str = "all",
    residue_type: str = "STY",
    radius: float = 6.0,
    residue_adjacent: int = 2,
    next_nearest: int = 3,
    ref_atom: str = "ca",
    file_extension: str = "h5",
    filename_format: str = "{to_process}-{residue_type}-{radius}A{residue_adjacent}R{next_nearest}N-{ref_atom}",

    verbose: bool = True,
) -> pd.DataFrame:
    """Load a motif dataframe from disk.


    Parameters
    ----------
    to_process : str, optional
        Which proteins to process, by default "all"
    residue_type : str, optional
        Which residues to process, by default "STY"
    radius : float, optional
        Radius of the motif, by default 6.0
    residue_adjacent : int, optional
        Number of sequence-adjacent residues the site
        of interest that are considered for
        - exclusion in the spatially-nearest neighbour 
            candidates
        - recording in the dataframe.
        Default value is `2`.
    next_nearest : int, optional
        Number of spatially-nearest neighbours to consider
        for each site of interest. Default value is `3`.
    ref_atom : str, optional
        Which atom to use for calculating spatial distances, 
        by default "ca" (alpha carbon)
    file_extension : str, optional
        File extension to use for loading the dataframe,
        by default "h5"
    filename_format : str, optional
        Format string for the filename, by default
        "{to_process}-{residue_type}-{radius}A{residue_adjacent}R{next_nearest}N-{ref_atom}"
    verbose : bool, optional
        Whether to print information about the loaded dataframe,    
        by default `True`.

    Returns
    -------
    pd.DataFrame
        Dataframe containing the motif data.

    """
    filename = filename_format.format(
        to_process=to_process,
        residue_type=residue_type,
        radius=int(radius),
        residue_adjacent=residue_adjacent,
        next_nearest=next_nearest,
        ref_atom=ref_atom,
    )
    outfile = DATA_DIR / "motif" / filename
    if file_extension == "h5":
        background_df = pd.read_hdf(outfile.with_suffix(f".{file_extension}"), key="df")
    elif file_extension == "csv":
        background_df = pd.read_csv(outfile.with_suffix(f".{file_extension}"))
    else:
        raise ValueError(f"Invalid file extension: {file_extension}")
    
    n_rows = len(background_df)
    n_proteins = len(background_df["protein_id"].unique())
    if verbose: print(f"Loaded {n_rows} rows from {n_proteins} proteins from {outfile}")
    return background_df

def format_motif_df(
    df: pd.DataFrame,
) -> pd.DataFrame:
    """Format a motif dataframe.

    Assumes 1-letter code for residues.
    
    """
    df["site_res"] = df["site"].apply(lambda x: x[0]) # Assumes 1-letter code
    return df


"""Exported variables."""
df = load_saved_motif_df()
df = format_motif_df(df)

# Filter for NaN 
df = df[~df["1_res"].isna()]
# all rows with +1 and -1 not nan
df = df[~df["-1"].isna()]
df = df[~df["+1"].isna()] 

background_df = df.copy()
del(df)

