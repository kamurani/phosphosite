
import numpy as np
import pandas as pd

from tqdm import tqdm
from pathlib import Path
from typing import Dict, List, Optional, Union, Callable

from phosphosite.structure.alphafold import make_af_dataframe, make_atomic_dataframe
from phosphosite.utils.structure import generate_sequence_from_df

from phosphosite.structure import StructureLoader
from phosphosite import AF_HUMAN_CIF
cif_loader = StructureLoader(AF_HUMAN_CIF, extension="cif.gz")

# TODO: join with phosphosite data. 

"""Process structural motif sites."""
def process(
    protein_ids: List[str],
    phosphosite_df: pd.DataFrame,
    reference_atom_dict: Dict[str, str] = {
        "S": "OG", 
        "T": "OG1",
        "Y": "OH",
    }, 
    residues_to_consider: str = "STY", 
    reference_atoms: List[str] = ["CA"],
    to_process: str = "all", # "p" for just phosphorylated.
    radius: float = 6.0,
    adjacency_range: int = 2,
    next_nearest: int = 3,
    verbose: bool = False,
    filepath: Optional[Union[str, Path]] = None,
    structure_loader: Callable = cif_loader,
) -> pd.DataFrame:
    """Process structural bubbles of residues around phosphosites.

    Parameters
    ----------
    protein_ids : List[str]
        List of protein ids to process.
    phosphosite_df : pd.DataFrame
        Phosphosite dataframe.
    reference_atom_dict : Dict[str, str], optional
        Dictionary mapping amino acid to reference atom, by default {"S": "OG", "T": "OG1", "Y": "OH"}
        The value is the cif code for the phosphorylated oxygen atom for the given residue.
    residues_to_consider : str, optional
        Residues to consider, by default "STY"
    reference_atoms : List[str], optional
        Reference atoms to consider in addition to those specified in `reference_atom_dict`,
        by default ["CA"]
    to_process : str, optional
        Whether to process all residues or just phosphorylated residues, by default "all". 
        One of "all" or "p".
    radius : float, optional
        Radius of bubble, by default 6.0
    adjacency_range : int, optional
        Range of adjacent residues to exclude from "next-nearest spatial
        neighbour" calculation, by default 2
        For example, if adjacent_range = 2, then residues -2 and +2 will be
        excluded from the candidates (relative to the site of interest's position).
    next_nearest : int, optional
        Number of next-nearest spatial neighbours to consider, by default 3.
        For example, if next_nearest = 3, then the 3 closest residues to the site
        of interest will be obtained.  If there is no residue within the radius, 
        then the value will be NaN.
    verbose : bool, optional
        Whether to print progress, by default False
    filepath : Optional[Union[str, Path]], optional
        Filepath to save dataframe, by default None
    structure_loader : Callable, optional
        Structure loader, by default `cif_loader` which loads .cif files from the 
        user-specified directory in the `config.yml` file.  
    
    Returns
    -------
    pd.DataFrame
        Dataframe of processed structural bubbles.
    """
    result_df = pd.DataFrame()

    for protein_id in protein_ids:
        atomic_df = make_atomic_dataframe(loader=cif_loader, protein_ids=[protein_id])
        sequence = generate_sequence_from_df(atomic_df)

        # Get list of all residues (position and amino acid) to consider. 
        # Foreach residue in residues_to_consider, get the position of the residue in the sequence.
        sites_to_process = [
            (res, pos) for res in residues_to_consider for pos in range(len(sequence)) if sequence[pos - 1] == res
        ]
        for res, pos in sites_to_process:
            try:
                row = atomic_df[(atomic_df["AA"] == res) & (atomic_df["position"] == pos)].iloc[0]   
            except IndexError:
                if verbose: tqdm.write(f"[{protein_id}] Could not find centre residue {res} at position {pos}")
                continue
            
            # Process bubble for this particular site. 
            records = []
            for ref_atom in reference_atoms + [reference_atom_dict[res]]:
                print(res, ref_atom)

                ref_coords = np.array(atomic_df[(atomic_df["AA"] == res) & (atomic_df["position"] == pos) & (atomic_df["atom_id"] == ref_atom)][["x_coord", "y_coord", "z_coord"]])
                if len(ref_coords) == 0:
                    if verbose: tqdm.write(f"[{protein_id}] No {ref_atom} atom for {res} at position {pos}")
                    print(reference_atom_dict)
                    continue
                
                # Filter by radius. 
                # Only include atoms (rows) whose distance from the reference atom is less than the radius.

                # Annotate with distance to reference atom.
                atomic_df[f"{ref_atom}_dist"] = np.linalg.norm(atomic_df[["x_coord", "y_coord", "z_coord"]].values - ref_coords, axis=1)

                # Filter to atoms within radius.
                dff = atomic_df[atomic_df[f"{ref_atom}_dist"] < radius]

                # Gamma Oxygen (and any ref atom)'s distance to itself should be zero.
                assert (dff[(dff.position == pos) & (dff.atom_id == ref_atom)][f"{ref_atom}_dist"] == 0).bool(), "yikes" 
                
                # Exclude residue positions. 
                to_exclude = list(range(pos - adjacency_range, pos + adjacency_range + 1))
                dff = dff[~dff["position"].isin(to_exclude)]

                # Get next-nearest spatial neighbours.
                nearest_neighbour = dff[dff[f"{ref_atom}_dist"] == dff[f"{ref_atom}_dist"].min()]

                if len(nearest_neighbour) == 0:
                    if verbose: tqdm.write(f"[{protein_id}] No nearest neighbour for {res} at position {pos}")
                    continue
                nearest_neighbour = nearest_neighbour.iloc[0]                

                prev_pos = pos - 1
                next_pos = pos + 1
                record = {
                    "protein_id": protein_id,
                    "prev": sequence[prev_pos - 1] if prev_pos > 0 else None, # 0-indexed
                    "res": res,
                    "pos": pos,
                    "next": sequence[next_pos - 1] if next_pos <= len(sequence) else None, # 0-indexed
                    "ref_atom": ref_atom,
                    "ref_atom_dist": nearest_neighbour[f"{ref_atom}_dist"],
                    "nn_res": nearest_neighbour["AA"],
                    "nn_pos": nearest_neighbour["position"],
                    "nn_atom": nearest_neighbour["atom_id"],
                }
                records.append(record)

            new_df = pd.DataFrame(records)
            result_df = pd.concat([result_df, new_df], axis=0)
             
    return result_df