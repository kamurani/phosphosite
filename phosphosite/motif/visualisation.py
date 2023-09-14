"""Visualisation of 3D motifs."""

from pathlib import Path 
import pandas as pd 
import numpy as np
import os 
import re 
import gzip 
import shutil
import Bio.PDB.MMCIF2Dict
from typing import Union, List, Tuple, Dict, Optional
from pathlib import Path

pd.options.mode.chained_assignment = None  # default='warn'

from phosphosite import DEFAULT_RADIUS
from phosphosite.utils import aa1to3, aa3to1
from phosphosite.utils.graphs import get_seq_distance

import plotly.express as px

def plot_heatmap(
    df: pd.DataFrame,
    aspect: str = None,
    title: str = f"Motif counts for nearest residues (R={DEFAULT_RADIUS}Å)",
    colour: str = "viridis", # "plasma" 
    height: int = 800,
    filepath: Path = None,
    show: bool = False,
    range_color: Tuple[float, float] = None,
):
    """Plot heatmap.
    
    Used for displaying motif counts. 
    """
    # df = df.T
    fig = px.imshow(
        df,
        color_continuous_scale=colour,
        # Don't enforce square 
        #width=1600,
        height=height,
        #labels=label_dict,
        title=title,
        aspect=aspect,
        range_color=range_color,
        
    )
    fig.update_xaxes(
        type='category', # In case we are using CIDs (integers)
        tickangle=45,
    ) 
    fig.update_yaxes(type='category')

    if filepath is not None:
        # Save png 
        if filepath.suffix == ".png":
            fig.write_image(str(filepath))
        # Save html
        elif filepath.suffix == ".html":
            fig.write_html(str(filepath))
        
    if show: fig.show()

    return fig    