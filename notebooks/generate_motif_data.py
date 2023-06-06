"""Script for generating motif data."""

import logging
logging.basicConfig(filename='example.log', encoding='utf-8', level=logging.WARNING)

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

from pathlib import Path 
import pandas as pd 
import numpy as np

from phosphosite.utils import aa1to3
from typing import Union, List, Tuple, Dict, Optional
from pathlib import Path

import graphein 
graphein.verbose(enabled=False)

N = 30000

def get_euc_dist(
    arr1: np.ndarray, arr2: np.ndarray
): 
    """Get euclidean distance between two arrays."""
    return np.sqrt(np.sum((arr1 - arr2) ** 2))
    

def get_node_id(
    site: str, 
    chain_id: str = "A",
) -> str: 
    mod_rsd, modification = site.split("-")
    aa = aa1to3[mod_rsd[0]]
    position = mod_rsd[1:]
    node_id = f"{chain_id}:{aa}:{position}"
    return node_id

def generate_node_id(
    node_dict: Dict[str, Union[str, int]],
    delimiter: str = ":",
) -> str: 
    return delimiter.join([str(node_dict[s]) for s in ["chain_id", "residue_name", "residue_number"]])


from phosphosite.dataset import phosphorylation # Filtered out isoforms
df = phosphorylation[phosphorylation["ORGANISM"] == "human"]

# Sort by ACC_ID
df = df.sort_values("ACC_ID")

# Filter by residue type, first character in MOD_RSD
allowed_residues = "STY"
df = df[df["MOD_RSD"].str[0].isin(list(allowed_residues))]

uniprot_id = "P51786"

af_version = 3
filename_template = "AF-{uniprot_id}-F1-model_v{af_version}.cif.gz"
filename = filename_template.format(uniprot_id=uniprot_id, af_version=af_version)

structure_dir = Path.home() / "STRUCTURAL_MOTIFS/DATA/"
af_cif_dir = structure_dir / "AF_HUMAN_CIF" 
af_pdb_dir = structure_dir / "AF_HUMAN_PDB"

# Assert that the cif and pdb directories exist.
assert af_cif_dir.exists()
assert af_pdb_dir.exists()


from phosphosite.structure import StructureLoader
cif_loader = StructureLoader(af_cif_dir, extension="cif.gz")
pdb_loader = StructureLoader(af_pdb_dir, extension="pdb")

protein_ids = list(df["ACC_ID"].unique())
existing_ids = cif_loader.get_existing_ids(protein_ids)
for protein_id in existing_ids:
    # Check that equivalent PDB file exists.
    if not pdb_loader.protein_id_exists(protein_id):
        print(f"Missing PDB file for {protein_id}")

df = df[df["ACC_ID"].isin(existing_ids)]


# Trim down uniprot ids for now. 
uniprot_ids = existing_ids[0:N]


"""Graphein subgraph method."""
pd.options.mode.chained_assignment = None  # default='warn'

from graphein.protein.config import ProteinGraphConfig
config = ProteinGraphConfig(verbose=False)

from graphein.protein.graphs import construct_graph
from graphein.protein.visualisation import plotly_protein_structure_graph


from phosphosite.graphs import get_motif_subgraph
from phosphosite.utils import aa1to3


failed = []
sequence_dict = {}

radius = 6.0
chain_id = "A"
force = False

out_dir = Path(f"data/processed/motifs")
out_file_format = "{uniprot_id}-R{radius}Å.csv"

skip_empty_motifs = False # if False, still write NaN to file.

for uniprot_id in uniprot_ids:

    # Check that .csv file exists for this uniprot_id.
    out_file = out_dir / out_file_format.format(
            uniprot_id=uniprot_id, 
            radius=radius,
        )
    if not force:
        if out_file.exists():
            print(f"Skipping {uniprot_id}. File exists.")
            continue

    # Construct graph. 
    pdb_path = pdb_loader.get_structure(uniprot_id)
    g = construct_graph(config=config, path=pdb_path)

    sequence_dict[uniprot_id] = g.graph[f"sequence_{chain_id}"]

    # Get sites for this uniprot_id
    site_df = df[df["ACC_ID"] == uniprot_id]
    sites = list(site_df["MOD_RSD"].unique())

    dict_list = []
    for site in sites:
        node_id = get_node_id(site)
        if node_id in g.nodes:
            s_g = get_motif_subgraph(
                g, 
                node_id, 
                radius = radius, # 5.0
            )
        else: 
            print(f"Node {node_id} not in graph ({uniprot_id}).")
            failed.append((uniprot_id, node_id))
            continue

        # Get nearest node to site.
        
        # Extract out adjacent residues (prev and next in sequence) 
        pos = int(node_id.split(":")[2])
        next_pos = pos + 1
        prev_pos = pos - 1
        seq = g.graph[f'sequence_{chain_id}']
        # sequence is 0 indexed 
        
        # Check if out of range. 
        if next_pos > len(seq):
            next_node_id = np.nan
        else:   
            next_node_id = f"{chain_id}:{aa1to3[seq[next_pos - 1]]}:{next_pos}"
        
        if prev_pos < 1:
            prev_node_id = np.nan
        else:
            prev_node_id = f"{chain_id}:{aa1to3[seq[prev_pos - 1]]}:{prev_pos}"

        site_coords = s_g.nodes[node_id]["coords"]

        candidate_nodes = [
            s_g.nodes[i]
            for i in s_g.nodes
            if i not in [prev_node_id, node_id, next_node_id] # Not the modified residue, or sequence adjacent residues.
        ]
        if len(candidate_nodes) == 0:
            print(f"No candidate nodes for {uniprot_id} {node_id}.")
            nearest_node = np.nan

            if skip_empty_motifs: 
                continue # don't create data for this motif.
        else:
            # Pick node with smallest euclidean distance to site.
            nearest_node = min(
                candidate_nodes,
                key=lambda node: get_euc_dist(site_coords, node["coords"]),
            )
            nearest_node = generate_node_id(nearest_node)
    
        dict_list.append(dict(
            uniprot_id=uniprot_id,
            prev=prev_node_id,
            site=node_id,
            next=next_node_id,
            nearest_node=nearest_node, # nearest that's prev or next. NaN if none. 
        ))
    
    
    motif_df = pd.DataFrame.from_dict(dict_list)
   
    motif_df.to_csv(
        out_file,
        index=False,
        # Tab delimited for easy reading.
        sep="\t",
        # Include column names.
        header=True,
    )