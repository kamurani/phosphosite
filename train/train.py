
"""Training script for phosphosite predictor."""

# TODO:
# - have 'setup' function that downloads the Pretrained embeddings from uniprot and puts it in correct directory etc. 
# - so we don't have to manually upload directory every time to cloud computing

import click as ck
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

# TODO
def generate_output_dataframe(data):
    return None



# TODO: 
# - run with 3 diff seeds to provide uncertainty estimates


@ck.command()
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
def main(
    dev: bool = False,
    model_name: str = "phosphogat",
    epochs: int = 200, 
):

    dropout = 0.1 
    batch_size = 64 # 100
    num_heads = 8 # 4 
    learning_rate = 0.001

    hidden_embedding_size = 256 

    device = "gpu" if torch.cuda.is_available() else "cpu"
    print(f"Using {device}.")

    learning_rate = 0.001

    

    np.random.seed(42)
    idx_all = np.arange(len(pyg_list))
    np.random.shuffle(idx_all)


    train_dataset 

    """Construct dataloaders."""
    from torch_geometric.data import DataLoader

    batch_size = 2 # for now
    data_list = 
    train_loader = DataLoader(data_list, batch_size=batch_size, shuffle=True, drop_last=True)
 



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
    

    trainer = Trainer(
        default_root_dir=SAVED_MODEL_DIR,
        max_epochs=epochs,
        gpus=1 if device == "gpu" else 0,
        callbacks=[early_stop_callback, checkpoint_callback],
        progress_bar_refresh_rate=1,
    )
    trainer.fit(model, train_loader, valid_loader)
    
    # evaluate on the model with the best validation set
    best_model = 

    


if __name__ == "__main__":
    main()