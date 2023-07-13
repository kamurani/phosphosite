import numpy as np 


def get_euc_dist(
    arr1: np.ndarray, arr2: np.ndarray
): 
    """Get euclidean distance between two arrays."""
    return np.sqrt(np.sum((arr1 - arr2) ** 2))