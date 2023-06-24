"""Processing structural data."""
import os 
import re 
import gzip 
import shutil
import Bio.PDB.MMCIF2Dict
import pandas as pd
import numpy as np 
import numba

from typing import Union, List, Tuple, Dict, Optional
from tqdm import tqdm 
from pathlib import Path

af_out_format = "AF-{uniprot_id}-F1-model_v{af_version}.{extension}"

"""Modified from structuremap 
https://github.com/MannLabs/structuremap 
"""
def process_af_data(
    directory: Union[str, Path], 
    protein_ids: Optional[List[str]] = None,
    out_format: str = None,
) -> pd.DataFrame:
    """Annotate alphafold structure data as a dataframe. 
    
    Modified from structuremap.processing.format_alphafold_data 
    [structuremap](https://github.com/MannLabs/structuremap)

    (modified to include .cif.gz files and custom file format)
    """
    alphafold_annotation_l = []
    protein_number = 0

    if out_format is None:
        out_format = "{uniprot_id}.cif"
    # Get list of uniprot_ids that we will use. 

    # NOTE: 
    # unsure how to 'extract' the uniprot_id from an arbitrary f-string. 
    # Usually you would use regex but the f string is just supplied as a string here. 
    # For now, we only support being passed uniprot_ids directly. (i.e. receiving 0 will not work)

    if protein_ids is None or len(protein_ids) == 0:
        # No protein_ids specified, so use all files in directory (that match out_format)
        pass 
    
    # Trim protein_ids down to only those that exist in the directory. 
    # NOTE: unnecessary if we check in the `for` loop.
    protein_ids = [
        protein_id 
        for protein_id in protein_ids
        if (directory / out_format.format(
            uniprot_id=protein_id,
        )).exists()
    ] 


    for protein_id in tqdm(protein_ids):

        filename = out_format.format(uniprot_id=protein_id) 
        filepath = directory / filename

        if not (filepath).exists(): continue 

        if filename.endswith("cif"):
            structure = Bio.PDB.MMCIF2Dict.MMCIF2Dict(filepath)
        elif filename.endswith("cif.gz"):
            with gzip.open(filepath, "rt") as f_in:
                structure = Bio.PDB.MMCIF2Dict.MMCIF2Dict(f_in)
        else: 
            raise ValueError(f"File '{filename}' has an unrecognized extension.")
            
        protein_number += 1
        df = pd.DataFrame({'protein_id': structure['_atom_site.pdbx_sifts_xref_db_acc'],
                            'protein_number': protein_number,
                            'AA': structure['_atom_site.pdbx_sifts_xref_db_res'],
                            'position': structure['_atom_site.label_seq_id'],
                            'quality': structure['_atom_site.B_iso_or_equiv'],
                            'atom_id': structure['_atom_site.label_atom_id'],
                            'x_coord': structure['_atom_site.Cartn_x'],
                            'y_coord': structure['_atom_site.Cartn_y'],
                            'z_coord': structure['_atom_site.Cartn_z']})

        df = df[df.atom_id.isin(['CA', 'CB', 'C', 'N'])].reset_index(drop=True)
        df = df.pivot(index=['protein_id',
                                'protein_number',
                                'AA', 'position',
                                'quality'],
                        columns="atom_id")
        df = pd.DataFrame(df.to_records())

        df = df.rename(columns={"('x_coord', 'CA')": "x_coord_ca",
                                "('y_coord', 'CA')": "y_coord_ca",
                                "('z_coord', 'CA')": "z_coord_ca",
                                "('x_coord', 'CB')": "x_coord_cb",
                                "('y_coord', 'CB')": "y_coord_cb",
                                "('z_coord', 'CB')": "z_coord_cb",
                                "('x_coord', 'C')": "x_coord_c",
                                "('y_coord', 'C')": "y_coord_c",
                                "('z_coord', 'C')": "z_coord_c",
                                "('x_coord', 'N')": "x_coord_n",
                                "('y_coord', 'N')": "y_coord_n",
                                "('z_coord', 'N')": "z_coord_n"})

        df = df.apply(pd.to_numeric, errors='ignore')
        df['secondary_structure'] = 'unstructured'
        if '_struct_conf.conf_type_id' in structure.keys():
            start_idx = [int(i) for i in structure['_struct_conf.beg_label_seq_id']]
            end_idx = [int(i) for i in structure['_struct_conf.end_label_seq_id']]
            note = structure['_struct_conf.conf_type_id']

            for i in np.arange(0, len(start_idx)):
                df['secondary_structure'] = np.where(
                    df['position'].between(
                        start_idx[i],
                        end_idx[i]),
                    note[i],
                    df['secondary_structure'])
        alphafold_annotation_l.append(df)

    alphafold_annotation = pd.concat(alphafold_annotation_l)
    alphafold_annotation = alphafold_annotation.sort_values(
        by=['protein_number', 'position']).reset_index(drop=True)

    alphafold_annotation['structure_group'] = [re.sub('_.*', '', i)
                                               for i in alphafold_annotation[
                                               'secondary_structure']]
    str_oh = pd.get_dummies(alphafold_annotation['structure_group'],
                            dtype='int64')
    alphafold_annotation = alphafold_annotation.join(str_oh)

    return(alphafold_annotation)


@numba.njit
def get_3d_dist(
    coordinate_array_1: np.ndarray,
    coordinate_array_2: np.ndarray,
    idx_1: int,
    idx_2: int
) -> float:
    """
    Function to get the distance between two coordinates in 3D space.
    Input are two coordinate arrays and two respective indices that specify
    for which points in the coordinate arrays the distance should be calculated.

    Parameters
    ----------
    coordinate_array_1 : np.ndarray
        Array of 3D coordinates.
        Must be 3d, e.g. np.float64[:,3]
    coordinate_array_2 : np.ndarray
        Array of 3D coordinates.
        Must be 3d, e.g. np.float64[:,3]
    idx_1 : int
        Integer to select an index in coordinate_array_1.
    idx_2 : int
        Integer to select an index in coordinate_array_2.

    Returns
    -------
    : float
        Distance between the two selected 3D coordinates.
    """
    dist = np.sqrt(
        (
            coordinate_array_1[idx_1, 0] - coordinate_array_2[idx_2, 0]
        )**2 + (
            coordinate_array_1[idx_1, 1] - coordinate_array_2[idx_2, 1]
        )**2 + (
            coordinate_array_1[idx_1, 2] - coordinate_array_2[idx_2, 2]
        )**2
    )
    return(dist)