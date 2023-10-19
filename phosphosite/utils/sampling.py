"""Functions for sampling a distribution of phosphosites."""

import pandas as pd 


from typing import Any, Union, List, Tuple, Dict, Optional

def sample_from_plddt_distribution(
    data: Any, 
    sample_from: Any, 
    secondary_structure_col: str, 
    # TODO:
    # - provide a list of discrete columns 
    # to use for bins 
    # in order to sample from a distribution across multiple columns
    plddt_col: str, 
    n_bins: int = 20,
    n_samples: int = None, 
    verbose: bool = False,

) -> pd.DataFrame: 
    """ 
    Generate a dataframe of negative examples by sampling from the distribution of pLDDT scores
    for a given (discrete) secondary structure annotation.

    Parameters
    ----------
    data: Any
        Dataframe containing the positive examples with the distribution to match. 
    sample_from: Any
        Dataframe containing the same `secondary_structure_col` and `plddt_col` columns as `data` but with
        a larger number of examples.  The distribution of pLDDT scores for each secondary structure annotation
        in `data` will be matched to the distribution of pLDDT scores for the same secondary structure annotation
        in `sample_from`.
    secondary_structure_col: str
        Name of the column containing the secondary structure annotation.
    plddt_col: str
        Name of the column containing the pLDDT score.
    n_bins: int
        Number of bins to use for the value counts.
    n_samples: int
        Number of samples to generate. If `None`, will generate the same number of samples as the positive examples.
        If a number greater than the number of positive examples is given, then the relative counts of the bins
        will be scaled up such that the total number of rows in the returned dataframe is equal to `n_samples`.

    Returns
    -------
    pd.DataFrame
        Dataframe containing the negative examples.
    """
    if n_samples is not None and n_samples <= len(sample_from):
        # Get ratio 
        ratio = n_samples / len(data)

    out_df = pd.DataFrame()
    for site_structure in data.site_structure.unique():
        pos_df = data[data[secondary_structure_col] == site_structure]
        pos_df = pos_df[pos_df[plddt_col].notna()]        

        candidates = sample_from[sample_from[secondary_structure_col] == site_structure]
        neg_df = candidates[candidates[plddt_col].notna()]    

        distr = pos_df[plddt_col].value_counts(bins=n_bins, sort=False)
        intervals = distr.index.to_list()
        bins = pd.IntervalIndex(intervals)

        # Assign each negative example to a bin (intervals)
        neg_df["binned"] = pd.cut(neg_df[plddt_col], bins)

        for bin in bins: 
            count = distr[bin] # value count from original positives

            # Scale count up to match in total 
            if n_samples is not None:
                count = int(count * ratio)
            if verbose: print(f"[{site_structure}] count: {count}")
            df = neg_df[neg_df.binned == bin]

            if len(df) == 0: # can't sample from empty dataframe with nonzero count 
                continue
            df = df.sample(n=count)
            out_df = out_df.append(df)

    return out_df
    


