import numpy as np 
import pandas as pd

def get_euc_dist(
    arr1: np.ndarray, arr2: np.ndarray
): 
    """Get euclidean distance between two arrays."""
    return np.sqrt(np.sum((arr1 - arr2) ** 2))


def generate_sequence_from_df(
    df: pd.DataFrame,
) -> str: 
    """Generate protein sequence from annotation dataframe.
    
    Assumes columns includes "position" and "AA". 
    """
    # Order by position. 
    dff = df.sort_values(by="position")["AA"].values
    dff = df.drop_duplicates(subset="position", keep="first")
    dff = dff.sort_values(by="position")["AA"].values 
    seq = "".join(dff)
    return seq