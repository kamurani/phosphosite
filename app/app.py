"""Dash app for phosphosite motif visualisation."""

from dash import Dash, dcc, html, Input, Output, State
import plotly.express as px

import pandas as pd
import numpy as np
from pathlib import Path

from phosphosite import GAMMA_OXYGEN_CODES
from phosphosite import DATA_DIR
from phosphosite import DEFAULT_RADIUS as radius
from phosphosite.utils import aa1to3, aa3to1
from phosphosite.motif.processing import make_count_df, make_motif_df
from phosphosite.motif.visualisation import plot_heatmap

from composition_stats import clr


# Load data
#from phosphosite.load_data import background_df as df
from phosphosite.bubble.data import result_df as df

app = Dash(__name__)

app.layout = html.Div([

    # Container
    html.Div([
        html.Div([
            
            # Control panel 
            html.Div([
                html.H1("Phosphosite motif visualisation"),
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
                        {"label": "Difference", "value": "diff"},
                        #{"label": "Entropy", "value": "entropy"},
                        {"label": "CLR", "value": "clr"},
                    
                    ],
                    value="none",
                    labelStyle={"display": "inline-block"},
                ),
                html.H3("Metric"),
                dcc.RadioItems(
                    id="scale", 
                    options=[
                        {"label": "Linear", "value": "linear"},
                        {"label": "Log2", "value": "log"},
                    ],
                    value="linear",
                    labelStyle={"display": "inline-block"},
                ),
                html.H3("Reference atom"),
                dcc.RadioItems(
                    id="ref_atom", 
                    options=[
                        {"label": "Alpha Carbon", "value": "CA"},
                        {"label": "Gamma Oxygen", "value": "oxygen"},
                    ],
                    value="CA",
                    labelStyle={"display": "inline-block"},
                ),
                html.H3("Group rows by"),
                dcc.RadioItems(
                    id="group_row_by", 
                    options=[
                        {"label": "Triplet", "value": "triplet"},
                        {"label": "-1", "value": "prev"},
                        {"label": "+1", "value": "next"},
                    ],
                    value="prev",
                    labelStyle={"display": "inline-block"},
                ),
                dcc.RadioItems(
                    id="normalise_3rd", 
                    options=[
                        {"label": "True", "value": True},
                        {"label": "False", "value": False},                
                    ],
                    value=True,
                    labelStyle={"display": "inline-block"},
                ),
            ], style={"display": "flex", "flex-direction": "column"}),

            # Distribution graph 
            html.Div([
                dcc.Graph(id="distribution"),
            ], style={"display": "flex"}),

        ], style={"display": "flex"}),

        # Two heatmaps side by side 
        html.Div([
            dcc.Graph(id="heatmap1"),
            dcc.Graph(id="heatmap2"),
            dcc.Graph(id="heatmap3"),
        ], style={"display": "flex"}),
    ]),
])

@app.callback(
    [
        Output("heatmap1", "figure"),
        Output("heatmap2", "figure"),
        Output("heatmap3", "figure"),
    ],
    [
        Input("residue", "value"),
        Input("normalisation", "value"),
        Input("scale", "value"),
        Input("ref_atom", "value"),
        Input("group_row_by", "value"),
        Input("normalise_3rd", "value")
    ],
)
def update_heatmaps(
    residue: str, 
    normalisation: str,
    scale: str,
    ref_atom: str,
    group_row_by: str,
    normalise_3rd: bool,

):
    """Update the heatmaps."""
    # Load data

    # Filter by residue 
    res_col = "site_res" if "site_res" in df.columns else "res"
    dff = df[df[res_col] == residue] 

    # Select reference atom
    if ref_atom == "CA":
        dff = dff[dff["ref_atom"] == ref_atom]
    elif ref_atom == "oxygen":
        dff = dff[dff["ref_atom"].isin(GAMMA_OXYGEN_CODES)]

    #dff.iloc["H-G", "R"] = 0

    psites = dff[dff["phosphosite"] == True]
    not_psites = dff[dff["phosphosite"] == False]

    version = "new"
    if version == "old":
        kwargs = dict(
            prev_col="-1",
            next_col="+1",
            nearest_col="1_res",
        )
    elif version == "new":
        kwargs = dict(
            prev_col="prev",
            next_col="next",
            nearest_col="nn_res",

            orient=group_row_by,
            #orient="prev",


        )


    placement = "vertical"
    if placement == "horizontal":
        motif_df = motif_df.T
        aspect = "auto"
        height = 1000
    else: 
        aspect = "auto"

        if group_row_by == "triplet": 
            height = 8000
        else:
            height = 1000
        



    plot_kwargs = dict(
        aspect=aspect,
        height=height,

    )
    
    def get_figure(
        dff: pd.DataFrame,
        title,
        **additional_kwargs,
    ):
        motif_df = make_motif_df(
            dff, 
            **kwargs,
            sep=residue,
        )
        fig = plot_heatmap(
            motif_df, 
            title=title,
            **plot_kwargs,
            **additional_kwargs,
            colour="plasma",
        )
        return fig

    
    # ALL SITES.
    fig1 = get_figure(
        dff, title=f"All {aa1to3[residue].capitalize()} residues in human proteome: next-nearest spatial neighbour (R={radius}Å)",
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
        psites, f"Phosphosite {aa1to3[residue].capitalize()} residues: next-nearest spatial neighbour (R={radius}Å)"
      
    )

    # NORMALISED
    motif_df_all = make_motif_df(
        dff, 
        **kwargs,
    )
    motif_df_p = make_motif_df(
        psites, 
        **kwargs,
    )

    # 
    

    # NORMALISED GRAPH.
    if normalisation is None or normalisation == "none": 
        normalised = motif_df_p / motif_df_all 
    elif normalisation == "all":
        motif_df_p = motif_df_p / motif_df_p.sum().sum()
        all_sites = motif_df_all / motif_df_all.sum().sum()
        normalised = motif_df_p / all_sites 
    elif normalisation == "row":
        # each row should sum to 1
        motif_df_p = motif_df_p.div(motif_df_p.sum(axis=1), axis=0)
        all_sites = motif_df_all.div(motif_df_all.sum(axis=1), axis=0)
        normalised = motif_df_p / all_sites

        
            # normalise by the third column. 
            # This is the - sum of (probability * log2(probability)) for each amino acid.
            # This is the entropy of each row. 
            # TODO
        
            

    elif normalisation == "diff":

        # subtract
        normalised = motif_df_p - motif_df_all 
        # log fold change per row.
    
    elif normalisation == "entropy":
        # Row normalisaton first 
        motif_df_p = motif_df_p.div(motif_df_p.sum(axis=1), axis=0)
        all_sites = motif_df_all.div(motif_df_all.sum(axis=1), axis=0)

        # Divide 
        normalised = motif_df_p / all_sites

        # Scale by amount of information in each row. 
        # This is the - sum of (probability * log2(probability)) for each amino acid.
        # This is the entropy of each row.
        normalised = normalised * motif_df_p.apply(lambda x: -np.sum(x * np.log2(x)), axis=1)
    elif normalisation == "clr":
        
        small_number = 1e-8
        # Add small value to avoid log(0)
        motif_df_p = motif_df_p + small_number
        motif_df_all = motif_df_all + small_number

        # Divide counts by counts 
        normalised = motif_df_p / motif_df_all

        # Apply clr 
        norm = clr(normalised.to_numpy()) 
        normalised = pd.DataFrame(norm, index=normalised.index, columns=normalised.columns)
        
    else: 
        raise ValueError(f"Unknown normalisation: {normalisation}")

    # ENTROPY.

    # Normalise the normalised heatmap 
    if normalise_3rd:
        # Make each row of the normalised matrix sum to 1. 
        normalised = normalised.div(normalised.sum(axis=1), axis=0) 

    if scale == "linear":
        pass
    elif scale == "log":
        normalised = np.log2(normalised)

    range_color = [0, 1]
    additional_kwargs = dict(
        #range_color=range_color,
    )
    fig3 = plot_heatmap(
        normalised, 
        title="Normalised frequency distributions of nearest spatial neighbour by triplet",
        **plot_kwargs,
        **additional_kwargs,
    )

    return fig1, fig2, fig3 

@app.callback(
    [
        Output("distribution", "figure"),
    ],
    [
        Input("residue", "value"),
        Input("ref_atom", "value"),
        Input("heatmap1", "hoverData"),
        Input("heatmap2", "hoverData"),
        Input("heatmap3", "hoverData"),
        Input("heatmap1", "clickData"),

    ]
)
def update_distribution_graph(
    residue: str,
    ref_atom: str,
    hover_1,
    hover_2,
    hover_3,
    click_1,

):
    print(f"click1: {click_1}")
    if hover_1 is not None:
        point = hover_1["points"][0]
        nn, motif = point["x"], point["y"]
    elif hover_2 is not None:
        point = hover_2["points"][0]
        nn, motif = point["x"], point["y"]
    else:
        nn, motif = "A", "ASX"

    print(f"hover 1", hover_1)
    print(f"hover2 ", hover_2)
    print(f"hover3", hover_3)

    # Filter dataframe 
    # Filter by residue 
    res_col = "site_res" if "site_res" in df.columns else "res"
    dff = df[df[res_col] == residue] 

    # Select reference atom
    if ref_atom == "CA":
        dff = dff[dff["ref_atom"] == ref_atom]
    elif ref_atom == "oxygen":
        dff = dff[dff["ref_atom"].isin(GAMMA_OXYGEN_CODES)]

    # CREATE DISTRIBUTION GRAPH 

    # "motif" column consists of 1st character of "prev", "res", "next" columns.
    dff["prev_res"] = dff["prev"].str[0]
    dff["next_res"] = dff["next"].str[0]

    if motif[0] == "X":
        dff["motif"] = "X" + dff["res"] + dff["next_res"]
    elif motif[2] == "X":
        dff["motif"] = dff["prev_res"] + dff["res"] + "X"
    else:
        dff["motif"] = dff["prev_res"] + dff["res"] + dff["next_res"]


    # Filter by motif and nearest neighbour 
    dff = dff[dff["motif"] == motif] 
    dff = dff[dff["nn_res"] == nn]

    

    # Create a histogram 
    
    distribution_graph = px.histogram(
        data_frame=dff,
        x="seq_dist",
        marginal="box",
        title=f"Distribution of sequence distance between {residue} sites and 3D-NN ({motif}-{nn}) (all sites)",
        nbins=1000,
        height=500,
        width=800,
    )
    
    
    
    return [distribution_graph]

if __name__ == "__main__":
    app.run_server(debug=True)

