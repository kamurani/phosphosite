"""AlphaFold 3D predicted structures."""
import Bio
import pandas as pd
import numpy as np
from pathlib import Path
from tqdm import tqdm
from typing import Union, List, Tuple, Dict, Optional

from phosphosite.structure.loader import StructureLoader
from phosphosite import AF_HUMAN_CIF as structure_dir


"""Create annotation dataframe of an alphafold structure."""
def make_af_dataframe(
    loader: StructureLoader,
    protein_ids: List[str],
) -> pd.DataFrame:
    """Annotate alphafold structure data as a dataframe. 
    
    Modified from structuremap.processing.format_alphafold_data 
    [structuremap](https://github.com/MannLabs/structuremap)

    (modified to include .cif.gz files and custom file format; 
    uses loader object instead of directory)
    """
    alphafold_annotation_l = []
    protein_number = 0

    # Get list of uniprot_ids that we will use. 

    if protein_ids is None or len(protein_ids) == 0:
        # No protein_ids specified, so use all files in directory (that match out_format)
        pass 
    
    # Trim protein_ids down to only those that exist in the directory. 
    # NOTE: unnecessary if we check in the `for` loop.
    protein_ids = [
        protein_id 
        for protein_id in protein_ids
        if loader.protein_id_exists(protein_id)
    ] 

    for protein_id in tqdm(protein_ids):

        structure = loader.parse_structure(protein_id)
            
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
