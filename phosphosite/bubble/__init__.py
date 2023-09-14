
import numpy as np
import pandas as pd

from tqdm.auto import tqdm
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
    result_df: Optional[pd.DataFrame] = None,
    overwrite: bool = False,
) -> pd.DataFrame:
    """Process structural bubbles of residues around phosphosites.

    Calculates next-nearest spatial neighbour residue based on reference 
    atom from a specified set of sites (e.g. phosphorylated residues).

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
    result_df : Optional[pd.DataFrame], optional
        Dataframe to append to, by default None
    overwrite : bool, optional
        Whether to overwrite existing file, by default False
    
    Returns
    -------
    pd.DataFrame
        Dataframe of processed structural bubbles.
    """
    if result_df is not None:
        tqdm.write("Using existing dataframe...")
        # Only process protein ids that have not been processed yet.
        protein_ids = list(set(protein_ids).difference(set(result_df["protein_id"].unique())))
        
        if verbose:
            tqdm.write(f"Loaded existing data from {filepath} with {len(result_df['protein_id'].unique())} uniprot ids.")
            if overwrite:
                tqdm.write("Overwriting existing data.")
            else:
                tqdm.write("Appending to existing data.")
            tqdm.write(f"Processing {len(protein_ids)} proteins.")
    else:
        result_df = pd.DataFrame()


    with tqdm(enumerate(protein_ids)) as pbar:
        for i, protein_id in pbar:
            pbar.update(1)
            pbar.set_description(f"{protein_id}")
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
                    
                    ref_coords = np.array(atomic_df[(atomic_df["AA"] == res) & (atomic_df["position"] == pos) & (atomic_df["atom_id"] == ref_atom)][["x_coord", "y_coord", "z_coord"]])
                    if len(ref_coords) == 0:
                        if verbose: tqdm.write(f"[{protein_id}] No {ref_atom} atom for {res} at position {pos}")
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
                        "site_qual": row["quality"],
                        "site_structure": row["secondary_structure"],
                        "nn_qual": nearest_neighbour["quality"],
                        "nn_structure": nearest_neighbour["secondary_structure"],
                    }
                    records.append(record)

                new_df = pd.DataFrame(records)
                result_df = pd.concat([result_df, new_df], axis=0)
                # Every 20 proteins, save to file.
                if i % 20 == 0:
                    save_dataframe(result_df, filepath)
            
            
    # Dump to file.
    save_dataframe(result_df, filepath)
    return result_df

def save_dataframe(df: pd.DataFrame, filepath: Path):
    if filepath is not None:
        if str(filepath).endswith(".h5"):
            df.to_hdf(filepath, key="data", mode="w")
        else:
            df.to_csv(filepath, index=False)


