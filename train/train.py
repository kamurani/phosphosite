
"""Training script for phosphosite predictor."""

# TODO:
# - have 'setup' function that downloads the Pretrained embeddings from uniprot and puts it in correct directory etc. 
# - so we don't have to manually upload directory every time to cloud computing

import click as ck
import glob
import re
import time
import numpy as np
import pandas as pd
import pytorch_lightning as pl
from pytorch_lightning.utilities.types import EVAL_DATALOADERS, STEP_OUTPUT
import torch as torch
import torch.nn as nn
import torch.nn.functional as F


from typing import Dict, List, Union
from pathlib import Path
from functools import reduce
from pytorch_lightning import Trainer
from torch.utils.data import DataLoader, Dataset
from pytorch_lightning.callbacks import StochasticWeightAveraging, EarlyStopping, ModelCheckpoint

from phosphosite import SAVED_MODEL_DIR
from phosphosite.model import PhosphoGAT

from phosphosite import PHOSPHOSITE_PREDICT_DIR
root_dir =  PHOSPHOSITE_PREDICT_DIR / "protein_graph_dataset"

# TODO
def generate_output_dataframe(data):
    return None



# TODO: 
# - run with 3 diff seeds to provide uncertainty estimates


@ck.command()
@ck.option(
    "--root-dir",
    "-r",
    # TODO: 
    # for now use str representation of path. 
    default=str(root_dir),
    help="Root directory for data.",
)
@ck.option(
    "--dev/--no-dev",
    default=False,
    help="Run a single training batch for debugging.",
)
@ck.option(
    "--model-name",
    "-n",
    default="M1",
    help="Name of the model to train.",
)
@ck.option(
    "--epochs",
    "-e",
    default=200,
    help="Number of epochs to train.",
)
@ck.option(
    "--dryrun/--no-dryrun",
    default=False,
    help="Run a dryrun of the training script.",
)
@ck.option(
    "--verbose/--no-verbose",
    default=True,
    help="Print more output.",
)
def main(
    root_dir: str = "",
    dev: bool = False,
    model_name: str = "phosphogat",
    epochs: int = 200, 
    dryrun: bool = False,
    verbose: bool = False,

):
    root_dir = Path(root_dir)
    if verbose: print(f"[Dataset] Using root directory: {root_dir}")
    dropout = 0.1 
    batch_size = 64 # 100
    num_heads = 8 # 4 
    learning_rate = 0.001

    hidden_embedding_size = 256 

    device = "gpu" if torch.cuda.is_available() else "cpu"
    print(f"Using {device}.")

    learning_rate = 0.001
    np.random.seed(42)
    #idx_all = np.arange(len(pyg_list))
    #np.random.shuffle(idx_all)

    """SETUP 
    - load graph construction functions etc.
    
    """
    from graphein.protein.config import ProteinGraphConfig
    from graphein.protein.edges.distance import add_distance_threshold
    from functools import partial

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
    from graphein.ml.conversion import GraphFormatConvertor

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

    """
    Create dataset.
    """
    # Load in the actual dataset (i.e. processed filenames)
    from phosphosite.ml.graph_dataset import PhosphositeGraphDataset 
    from phosphosite import INDEX_DICT_PATH
    
    processed_filenames = [Path(a).stem for a in glob.glob(str(root_dir / "processed" / "*.pt"))]
    if verbose: print(f"Using {len(processed_filenames)} processed files.")


    from phosphosite.utils.io import load_index_dict
    indexes_dict = load_index_dict(filepath=INDEX_DICT_PATH)

    kwargs = dict(
        root=root_dir,
        graphein_config=config, 
        graph_transformation_funcs=graph_transforms,
        graph_format_convertor=convertor,
        pre_transform=None, # before saved to disk , after PyG conversion 
        pre_filter=None,    # whether it will be in final dataset
    )
    ds = PhosphositeGraphDataset(
        uniprot_ids=processed_filenames,
        y_label_map=indexes_dict,
        **kwargs,
    )
    if verbose: print(ds)


    from phosphosite.ml import get_dataloader_split
    train_loader, valid_loader, test_loader = get_dataloader_split(ds, batch_size=32, train_batch_size=32)
    if verbose: print(train_loader, valid_loader, test_loader)

 
    model = PhosphoGAT(
        dropout=dropout,
        batch_size=batch_size,
        learning_rate=learning_rate,
    )
    """Train model."""
    # Early stopping 
    patience = 10
    early_stop_callback = EarlyStopping(
        monitor="val_loss",
        min_delta=0.00,
        patience=patience,   
        verbose=True,
        mode="min",
    )

    checkpoint_callback = ModelCheckpoint(
        monitor="val_loss",
        dirpath=SAVED_MODEL_DIR,
        filename="model-{epoch:02d}-{val_loss:.2f}",
        save_top_k=1,
        mode="min",
    )
    

    trainer = Trainer(
        default_root_dir=SAVED_MODEL_DIR,
        max_epochs=epochs,
        accelerator="auto", devices="auto",
        callbacks=[early_stop_callback, checkpoint_callback],
        enable_progress_bar=True,
        fast_dev_run=dev,
    )
    trainer.fit(model, train_loader, valid_loader)
    
    # evaluate on the model with the best validation set
    trainer.test(ckpt_path="best", test_dataloaders=test_loader)

    


if __name__ == "__main__":
    main()