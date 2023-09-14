"""ML models for predicting phosphorylation sites."""

# pytorch geometric and pytorch lightning
from typing import Any
import torch
import torch.nn as nn
import torch.nn.functional as F

from torch.utils.data import Dataset, DataLoader

from torch_geometric.nn import GCNConv, GATConv, GATv2Conv, global_add_pool, global_mean_pool, global_max_pool, Sequential as PyGSequential

# pytorch lightning
import pytorch_lightning as pl

# import metrics from torch
#from torchmetrics import Accuracy, Precision, Recall, F1, AUROC, ConfusionMatrix


from phosphosite.model.defaults import (
    SEQUENCE_EMBEDDING_SIZE,
)



class PhosphoGAT(pl.LightningModule):
    """Graph Attention Network for predicting phosphorylation sites.
    
    Parameters
    ----------

    """
    def __init__(
        self, 

        # Model params
        node_embedding_size: int = SEQUENCE_EMBEDDING_SIZE, 
        num_heads: int = 8, 
        dropout: float = 0.1,

        # Training params 
        learning_rate: float = 0.001,
        batch_size: int = 32,
        num_workers: int = 8,
        num_epochs: int = 200,


        train_dataset: Dataset = None,
        validation_dataset: Dataset = None,
        test_dataset: Dataset = None,

    ) -> None:
        super().__init__()
        self.learning_rate = learning_rate
        self.batch_size = batch_size
        self.num_workers = num_workers
        self.num_epochs = num_epochs

        self.train_dataset = train_dataset
        self.validation_dataset = validation_dataset
        self.test_dataset = test_dataset

        self.out_heads = 1 

        hidden_size = node_embedding_size // num_heads  

        out_size = hidden_size*num_heads # change this to be less perhaps?

        # TODO: 
        # - test versions of allowing more or less dimensions (128 up to 1024 (same as original embeddings)) 
        # in the output layers of the Conv

        """Layers"""
        self.conv1 = GATv2Conv(node_embedding_size, hidden_size, heads=num_heads, dropout=dropout)
        self.conv2 = GATv2Conv(hidden_size*num_heads, out_size, heads=self.out_heads, dropout=dropout)

        self.convolutions = PyGSequential(
            "x, edge_index", [
                # Dropout before the conv layers? Hmm.
                (nn.Dropout(dropout), "x -> x"),
                (self.conv1, "x, edge_index -> x"),
                nn.SELU(),
                nn.Dropout(dropout),
                (self.conv2, "x, edge_index -> x"),
                nn.SELU(),
            ]
        )
        # After the conv layers, we should ideally have a learnt representation of the structural / sequence 
        # information of the protein per each site. We can then apply the fully connected layer to each
        # node to get a phosphosite prediction. 

        # TODO: 
        # use b_factor as a kind of mask? 

        """Classification"""

        # Several fully connected layers with dropout and selu 
        self.classifier = nn.Sequential(
            nn.Linear(out_size, 512),
            nn.SELU(),
            nn.Dropout(dropout),
            nn.Linear(512, 256),
            nn.SELU(),
            nn.Dropout(dropout),
            nn.Linear(256, 128),
            nn.SELU(),
            nn.Dropout(dropout),
            nn.Linear(128, 1),
            # binary classification i.e. sigmoid
            nn.Sigmoid(),
        )

    def forward(
        self, 
        protein_data, 

    ) -> torch.Tensor:
        """Forward pass of the model.

        Parameters
        ----------
        protein_data : torch_geometric.data.Data
            The input data object.

        Returns
        -------
        torch.Tensor
            The output tensor of shape (batch_size, output_dim)

        """
        
        x, edge_index, batch = protein_data.x, protein_data.edge_index, protein_data.batch


        x = self.convolutions(x, edge_index)
        #pooled = global_add_pool(x, batch)
       
        # Apply fully connected layers to each node
        x = self.classifier(x)
        return x



