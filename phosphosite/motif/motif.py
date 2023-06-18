"""Structural motif data."""

from pathlib import Path 
import pandas as pd 
import numpy as np

# Import structuremap functions
import structuremap.utils
structuremap.utils.set_logger()
from structuremap.processing import download_alphafold_cif, download_alphafold_pae, format_alphafold_data, annotate_accessibility, get_smooth_score, annotate_proteins_with_idr_pattern, get_extended_flexible_pattern, get_proximity_pvals, perform_enrichment_analysis, perform_enrichment_analysis_per_protein, evaluate_ptm_colocalization, extract_motifs_in_proteome
from structuremap.plotting import plot_enrichment, plot_ptm_colocalization

import os 
import re 
import gzip 
import shutil
import Bio.PDB.MMCIF2Dict
from typing import Union, List, Tuple, Dict, Optional
from pathlib import Path

from phosphosite.structure.processing import process_af_data

from phosphosite import DATA_DIR
annotation_dir = DATA_DIR / "structure_annotations"


"""Initialise phosphosite dataset."""
from phosphosite.dataset import phosphorylation # Filtered out isoforms
df = phosphorylation[phosphorylation["ORGANISM"] == "human"]

# Sort by ACC_ID
df = df.sort_values("ACC_ID")

# Filter by residue type, first character in MOD_RSD
allowed_residues = "STY"
df = df[df["MOD_RSD"].str[0].isin(list(allowed_residues))]

"""Initialise structure loader."""
af_version = 3
filename_template = "AF-{uniprot_id}-F1-model_v{af_version}.cif.gz"

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

# only use existing 
df = df[df["ACC_ID"].isin(existing_ids)]

N = 10000
from phosphosite.structure.processing import process_af_data
structure_df = process_af_data(
    af_cif_dir, 
    out_format="AF-{uniprot_id}-F1-model_v3.cif.gz",
    protein_ids=existing_ids, 
    #protein_ids=existing_ids[0:N],
)

# Save 
filepath = annotation_dir / f"structure_df_{N}.csv"
structure_df.to_csv(filepath, sep="\t", index=False)



exit()
if filepath.exists():
    structure_df = pd.read_csv(filepath, sep="\t")
    todo_ids = [x for x in uniprot_ids if x not in structure_df["protein_id"].unique()]
    # New batch 
    batch_df = process_af_data(
        af_cif_dir, 
        out_format="AF-{uniprot_id}-F1-model_v3.cif.gz",
        #protein_ids=existing_ids, 
        protein_ids=todo_ids,
    )
    print(f"length of batch_df: {len(batch_df)}")
    # Concat batch and existing 
    structure_df = pd.concat([structure_df, batch_df], axis=0)
    # Save 
    structure_df.to_csv(filepath, sep="\t", index=False)

else: 
    print(f"Initialising structure_df")
    structure_df = process_af_data(
        af_cif_dir, 
        out_format="AF-{uniprot_id}-F1-model_v3.cif.gz",
        #protein_ids=existing_ids, 
        protein_ids=uniprot_ids[0:10],
    )

    # Save 
    structure_df.to_csv(filepath, sep="\t", index=False)




