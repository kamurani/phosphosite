"""Working with features on graphs."""

import networkx as nx
import numpy as np


from typing import List, Dict, Tuple, Union, Any

def add_residue_embedding(
    G: nx.Graph,
    embedding: np.ndarray,
    label: str = "embedding",
) -> nx.Graph:
    """
    Adds residue embedding to graph.

    Assumes that the embedding provided has the same 
    dimension 1 as the number of nodes in the graph.

    NOTE: works with one chain only.
    
    """
    if embedding.shape[0] != len(G.nodes):
        raise ValueError(f"Embedding shape {embedding.shape} does not match graph size {len(G.nodes)}.")

    # TODO: chain selection?

    for i, (n, d) in enumerate(G.nodes(data=True)):
        G.nodes[n][label] = embedding[i]

    return G