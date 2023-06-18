"""Processing dataframes for motif analysis."""


from pathlib import Path 
import pandas as pd 
import numpy as np
from typing import Union, List, Tuple, Dict, Optional
from pathlib import Path

pd.options.mode.chained_assignment = None  # default='warn'

from phosphosite import DATA_DIR
from phosphosite.utils import aa1to3, aa3to1



def make_count_df(
    df: pd.DataFrame,
):
    count_df = df.groupby(["prev", "next"]).size().reset_index(name="count")
    # Turn into a 2d matrix where each row is a prev and each column is a next.
    count_df = count_df.pivot_table(
        index="prev",
        columns="next",
        values="count",
        fill_value=0,
    )
    return count_df

def make_motif_df(
    df: pd.DataFrame,
    prev_col: str = "prev",
    next_col: str = "next",
    nearest_col: str = "nearest_node",
) -> pd.DataFrame:
    """Return a dataframe with the counts of each residue for each motif.

    Parameters
    ----------
    df : pd.DataFrame
        Dataframe containing the columns prev_col, next_col, nearest_col.
    prev_col : str, optional
        Column name for the previous residue, by default "prev"
    next_col : str, optional
        Column name for the next residue, by default "next"
    nearest_col : str, optional 
        Column name for the nearest residue, by default "nearest_node"

    Returns
    -------
    pd.DataFrame
        Dataframe with the counts of each residue for each motif.
    """
    def get_residue(node: str) -> str:
        """Return the residue from a node."""

        if ":" in node:
            # A:SER:123 FORMAT
            return node.split(":")[1]
        else:
            # S123 FORMAT 
            return node[0]

    df["nearest_res"] = df[nearest_col].apply(lambda x: get_residue(x))

    # Turn prev and next into just residue. 

    for col in [prev_col, next_col]:
        df[col] = df[col].apply(lambda x: get_residue(x) if x is not np.nan else x)
    # For each combination of (prev, next); count the number of times each nearest_res occurs. 
    motif_counts = df.groupby([prev_col, next_col, "nearest_res"]).size().reset_index(name="count")
    # Collapse to one row per combination of (prev, next), and a column for each nearest_res.  
    # If NaN, replace with 0.
    
    motif_counts = motif_counts.pivot_table(
        index=[prev_col, next_col],
        columns="nearest_res",
        values="count",
        fill_value=0,
    ).reset_index()
    motif_counts
    # Make the prev and next columns into a single column and use as the index.
    motif_counts["motif"] = motif_counts.apply(
        lambda row: f"{row[prev_col]}-{row[next_col]}",
        axis=1,
    )
    motif_counts = motif_counts.set_index("motif")

    # Drop the prev and next columns.
    motif_counts = motif_counts.drop(columns=[prev_col, next_col])
    return motif_counts


def get_top_hit_df(motif_df: pd.DataFrame) -> pd.DataFrame:
    """Return a dataframe with the top hit for each motif."""
    # New column "res" which stores the column name of the largest value in each row.
    df = pd.DataFrame()
    df.index = motif_df.index
    df["res"] = motif_df.idxmax(axis=1)
    
    # New column "count" which stores the largest value in each row.
    df["count"] = motif_df.max(axis=1)
    return df

def make_adjacency_matrix(
    motif_df: pd.DataFrame,
    column: str = "count",
) -> pd.DataFrame:
    """Turn df into an adjacency matrix.
    
    Assumes index is of form "A-B"
    """
    # Split the index into prev and next columns.
    motif_df["prev"] = motif_df.index.str.split("-").str[0]
    motif_df["next"] = motif_df.index.str.split("-").str[1]
    # Turn into a 2d matrix where each row is a prev and each column is a next.


    # Check if column is numeric. 
    if motif_df[column].dtype == "object":

        #raise ValueError(f"Column {column} is not numeric.")
        
        # Make a dictionary of the unique values in the column.
        unique_values = motif_df[column].unique()
        value_dict = {value: i for i, value in enumerate(unique_values)}
        # Replace the values in the column with the index of the value in the dictionary.
        motif_df[column] = motif_df[column].apply(lambda x: value_dict[x])
        df = motif_df.pivot_table(
            index="prev",
            columns="next",
            values=column,
            fill_value=0,
        )
        
        # Convert back to the original values.
        df = df.applymap(lambda x: unique_values[x])
        return df
    else:
        df = motif_df.pivot_table(
            index="prev",
            columns="next",
            values=column,
            fill_value=0,
        )
        return df
