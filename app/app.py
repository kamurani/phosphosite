"""Dash app for phosphosite motif visualisation."""

from dash import Dash, dcc, html, Input, Output, State
import plotly.express as px

import pandas as pd
import numpy as np
from pathlib import Path

from phosphosite import DATA_DIR
from phosphosite import DEFAULT_RADIUS as radius
from phosphosite.utils import aa1to3, aa3to1
from phosphosite.motif.processing import make_count_df, make_motif_df
from phosphosite.motif.visualisation import plot_heatmap

# Load data
from phosphosite.load_data import background_df as df

app = Dash(__name__)

app.layout = html.Div([
    html.H1("Phosphosite motif visualisation"),
    html.Div([
        html.H3("Select residue"),
        
        # Which residue to select. Radio button 
        dcc.RadioItems(
            id="residue", 
            options=[
                {"label": "Serine", "value": "S"},
                {"label": "Threonine", "value": "T"},
                {"label": "Tyrosine", "value": "Y"},
            ],
            value="S",
            labelStyle={"display": "inline-block"},
        ),
        html.H3("Normalisation"),
        dcc.RadioItems(
            id="normalisation", 
            options=[
                {"label": "None", "value": "none"},
                {"label": "All", "value": "all"},
                {"label": "Row", "value": "row"},
                {"label": "Log fold change per row", "value": "diff"},
              
            ],
            value="none",
            labelStyle={"display": "inline-block"},
        ),

        # Two heatmaps side by side 
        html.Div([
            dcc.Graph(id="heatmap1"),
            dcc.Graph(id="heatmap2"),
        ], style={"display": "flex"}),
    ]),
])

@app.callback(
    [
        Output("heatmap1", "figure"),
        Output("heatmap2", "figure"),
    ],
    [
        Input("residue", "value"),
        Input("normalisation", "value"),
    ],
)
def update_heatmaps(
    residue: str, 
    normalisation: str,
):
    """Update the heatmaps."""
    # Load data

    # Filter by residue
    dff = df[df["site_res"] == residue] 
    psites = dff[dff["phosphosite"] == True]
    not_psites = dff[dff["phosphosite"] == False]

    kwargs = dict(
        prev_col="-1",
        next_col="+1",
        nearest_col="1_res",
    )

    placement = "vertical"
    if placement == "horizontal":
        motif_df = motif_df.T
        aspect = "auto"
        height = 1000
    else: 
        aspect = "auto" 
        height = 8000


    plot_kwargs = dict(
        aspect=aspect,
        height=height,
    )
    
    def get_figure(
        dff: pd.DataFrame,
        title,
    ):
        motif_df = make_motif_df(
            dff, 
            **kwargs,
        )
        fig = plot_heatmap(
            motif_df, 
            title=title,
            **plot_kwargs,
        )
        return fig
    
    # ALL SITES.
    fig1 = get_figure(
        dff, title=f"ALL {aa1to3[residue].capitalize()} residues in human proteome: next-nearest spatial neighbour (R={radius}Å)",
    )

    # PHOSPHOSITES.
    if False:
        if normalisation is None or normalisation == "none": 
            pass 
        elif normalisation == "all":
            psites = psites / psites.sum().sum()
        elif normalisation == "row":
            # each row should sum to 1
            psites = psites.div(psites.sum(axis=1), axis=0)
        elif normalisation == "diff":
            # first, normalise background. 
            dff = dff / dff.sum().sum()

            # turn psites into frequencies. 
            psites = psites / psites.sum().sum()

            # log fold change per row.
            psites = np.log2(psites / dff)
        else: 
            raise ValueError(f"Unknown normalisation: {normalisation}")

    # first, normalise background. 
    #dff = dff / dff.sum().sum()

    # turn psites into frequencies. 
    #psites = psites / psites.sum().sum()

    # log fold change per row.
    #psites = np.log2(psites / dff)

    fig2 = get_figure(
        psites, f"PHOSPHOSITE {aa1to3[residue].capitalize()} residues: next-nearest spatial neighbour (R={radius}Å)"
      
    )


    return fig1, fig2


if __name__ == "__main__":
    app.run_server(debug=True)

