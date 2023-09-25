"""Generate PyG graph."""



import torch_geometric
import networkx as nx 

from pathlib import Path
from functools import partial
from typing import Dict, List, Union, Any, Optional, Tuple, Callable
from phosphosite.graphs import (
    NODE_DISTANCE_THRESHOLD,
    LONG_INTERACTION_THRESHOLD,

)

from graphein.protein.config import ProteinGraphConfig
from graphein.protein.graphs import construct_graph

from graphein.protein.edges.distance import add_distance_threshold

from graphein.ml.conversion import GraphFormatConvertor

from phosphosite.protein.embeddings import get_embedding 
from phosphosite.graphs.features import add_residue_embedding

def add_per_residue_embedding(
    g: nx.Graph, 
    verbose: bool = False,
    uniprot_id_attribute: str = "name",
) -> nx.Graph:
    """
    Retrieve per-residue embedding from graph.

    Uses uniprot_id to retrieve embedding from database.

    Parameters
    ----------
    g : nx.Graph
        Graph to add embedding to.
    verbose : bool, optional
        Whether to print verbose output, by default False
    uniprot_id_attribute : str, optional
        Name of the attribute in the original graph that contains 
        the uniprot_id, by default "name".
    
    Returns
    -------
    nx.Graph
        Graph with per-residue embedding added.
    
    Raises
    ------
    ValueError
        If uniprot_id is not in graph.


    """
    uniprot_id = g.__getattribute__(uniprot_id_attribute)
    # add embedding
    try:
        emb = get_embedding(uniprot_id)
    except ValueError as e:
        if verbose: print(f"{uniprot_id}: {e}")
        return None

    try:
        g = add_residue_embedding(g, emb, label="x")
    except ValueError as e:
        if verbose: print(f"{uniprot_id}: {e}")
        return None
    
    return g

def get_pyg_graph(
    uniprot_id: str,
    pdb_dir: Path,
    filename_format: str = "{uniprot_id}.pdb", #"AF-{uniprot_id}-F1-model_v4.pdb", 
    edge_threshold_distance: float = NODE_DISTANCE_THRESHOLD,
    long_interaction_threshold: int = LONG_INTERACTION_THRESHOLD,


    **kwargs,
) -> torch_geometric.data.Data:
    """Get PyG graph for given uniprot_id.

    Uses threshold 

    Parameters
    ----------
    uniprot_id : str
        Uniprot ID of protein.
    pdb_dir : Path
        Path to directory containing PDB files.
    filename_format : str, optional
        Format of PDB filename, by default "AF-{uniprot_id}-F1-model_v4.pdb"
    edge_threshold_distance : float, optional
        Threshold distance for edges, by default NODE_DISTANCE_THRESHOLD
    long_interaction_threshold : int, optional
        Threshold for long-range interactions, by default LONG_INTERACTION_THRESHOLD

    Returns
    -------
    nx.Graph
        Graph of protein structure.

    """

    # Setup config
    new_edge_funcs = {"edge_construction_functions": [
        partial(
        add_distance_threshold, long_interaction_threshold=long_interaction_threshold, threshold=edge_threshold_distance)
    ]}
    config = ProteinGraphConfig(
        pdb_dir=pdb_dir,

        granularity="CA",

        # Node features
        #node_metadata_functions=[],

        # Edges based on thresholded distance 
        **new_edge_funcs,
    )

    """Create graph"""
    pdb_path = Path(pdb_dir) / filename_format.format(uniprot_id=uniprot_id)
    g = construct_graph(config=config, path=pdb_path, verbose=False)

    # add embedding
    try:
        emb = get_embedding(uniprot_id)
    except ValueError as e:
        print(f"{uniprot_id}: {e}")
        return None

    try:
        g = add_residue_embedding(g, emb, label="x")
    except ValueError as e:
        print(f"{uniprot_id}: {e}")
        return None

    
    columns = [
        "b_factor",
        #"coords",
        "name",
        "edge_index",
        "x", # T5 per-residue embedding

    ]
    convertor = GraphFormatConvertor(
        src_format="nx", dst_format="pyg", verbose="gnn",
        columns=columns,
    )
    pyg = convertor(g)
    assert type(pyg) is torch_geometric.data.Data
    assert pyg.validate(raise_on_error=True) == True, "PyG graph is invalid."

    return pyg

