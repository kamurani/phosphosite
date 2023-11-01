import torch
import glob

from typing import Dict, List, Union, Any

from phosphosite.ml.graph_dataset import PhosphositeGraphDataset 
from pathlib import Path
from tqdm import tqdm
from functools import partial
from graphein.protein.config import ProteinGraphConfig
from graphein.protein.edges.distance import add_distance_threshold
from graphein.ml.conversion import GraphFormatConvertor






from phosphosite import PHOSPHOSITE_PREDICT_DIR 

def get_dataset(
    root_dir: Path = PHOSPHOSITE_PREDICT_DIR / "protein_graph_dataset",
    index_dict: Any = None,
    uniprot_ids: List[str] = None, 
):
    if uniprot_ids is None: 
        # Use just processed files in existence.
        processed_filenames = [Path(a).stem for a in glob.glob(str(root_dir / "processed" / "*.pt"))]
        uniprot_ids = processed_filenames

    long_interaction_threshold = 5 # seq positions 
    edge_threshold_distance = 6.0 # Å
    new_edge_funcs = {"edge_construction_functions": [
        partial(
        add_distance_threshold, long_interaction_threshold=long_interaction_threshold, threshold=edge_threshold_distance)
    ]}
    config = ProteinGraphConfig(
        granularity="CA",
        **new_edge_funcs,
    )

    columns = [
        "b_factor",
        "name",
        "edge_index",
        "x", # T5 per-residue embedding
    ]
    convertor = GraphFormatConvertor(
        src_format="nx", dst_format="pyg", verbose="gnn",
        columns=columns,
    )

    # List of functions that consume a nx.Graph and return a nx.Graph. Applied to graphs after construction but before conversion to pyg
    from phosphosite.graphs.pyg import add_per_residue_embedding
    graph_transforms = [
        add_per_residue_embedding,
    ]

    kwargs = dict(
        root=root_dir,
        graphein_config=config, 
        graph_transformation_funcs=graph_transforms,
        graph_format_convertor=convertor,
        pre_transform=None, # before saved to disk , after PyG conversion 
        pre_filter=None,    # whether it will be in final dataset
    )
    ds = PhosphositeGraphDataset(
        uniprot_ids=uniprot_ids,
        y_label_map=index_dict,
        **kwargs,
    )
    return ds

dataset = get_dataset()