import requests 
import time 
import pandas as pd 
import numpy as np
import Bio

from typing import List
from pathlib import Path
from Bio import SwissProt
from Bio import ExPASy

from phosphosite import UNIPROT_DATA_DIR




# Using Bio.SwissProt, download the record for the protein of interest.
#



def download_uniprot_data(
    protein_ids: List[str],
    data_dir: Path = UNIPROT_DATA_DIR,
    url: str = "https://rest.uniprot.org/uniprotkb/{protein_id}.txt",
) -> None:
    """Retrieve entries from uniprot KB"""
    for protein in protein_ids:
        download_url = url.format(protein_id=protein)

        # Save to file in data directory
        filepath = data_dir / f"{protein}.txt"
        if not filepath.exists():
            r = requests.get(download_url)
            with open(filepath, "wb") as f:
                f.write(r.content)

        # Sleep to avoid overloading server
        time.sleep(0.5)
            



"""Query realtime to get record"""
def get_record(protein_id: str) -> SwissProt.Record:
    """Get record from UniProt for a given protein id."""
    handle = ExPASy.get_sprot_raw(protein_id)
    record = SwissProt.read(handle)
    return record


def map_site_to_domain(
    protein_id: str, 
    res: str = None,
    pos: int = None,
    uniprot_data_dir: Path = UNIPROT_DATA_DIR,
    file_name_format: str = "{protein_id}.txt",
    output_format: str = "df", # None
    multiple: bool = False,
) -> str:
    """Given a protein id, residue, and position, return the domain that the site is in.
    
    Uses the uniprot KB text download into the uniprot_data_dir, and downloads it if it doesn't exist.
    
    Parameters
    ----------
    protein_id : str
        Uniprot protein id.
    res : str
        Residue (single letter code).
    pos : int
        Position of residue.
    uniprot_data_dir : Path, optional
        Path to uniprot data directory, by default UNIPROT_DATA_DIR
    
    Returns
    -------
    str
        Domain that the site is in.
    """
    filename = file_name_format.format(protein_id=protein_id)
    filepath = uniprot_data_dir / filename
    if not filepath.exists():
        download_uniprot_data(protein_ids=[protein_id], data_dir=uniprot_data_dir)
    
    with open(filepath, "r") as f:
        r = SwissProt.parse(f)
    
        # Yield the record.
        r = next(r)
        
        row = {}
        for f in r.features:
            # If unknown, skip.

            if type(f.location.start) == Bio.SeqFeature.UnknownPosition or type(f.location.end) == Bio.SeqFeature.UnknownPosition:
                # Avoid error when comparing to pos. 
                continue
            
            if pos in f.location:
                
                
                if f.type == "DOMAIN":
                    row["domain"] = f.qualifiers["note"]
                

                for feature in ["ZN_FING", "COILED", "DNA_BIND"]:
                    if f.type == feature:
                        row[feature] = True
                    continue

            
       
        
        return row
                
from phosphosite import GAMMA_OXYGEN_CODES

def get_proteins_from_motif(
    df: pd.DataFrame,
    triplet: str, 
    nn: str, 
    ref_atom: str = None,
    phosphosite: bool = None,
) -> List[str]:
    """Return the list of proteins that contain a given motif."""
    
    # Filter by ref_atom
    if ref_atom is not None:
        if ref_atom == "CA": ref_atom = ["CA"]
        elif ref_atom == "oxygen": ref_atom = GAMMA_OXYGEN_CODES
        df = df[df["ref_atom"].isin(ref_atom)]
    
    # Turn dataframe cols into just residue (i.e. first character)
    for col in ["prev", "res", "next", "nn_res"]:
        df[col] = df[col].str[0] 

    # If "next" is unspecified: 
    #if len(triplet) == 2:
    #    df = df[(df["prev"] == triplet[0]) & (df["res"] == triplet[1]) & (df["nn_pos"] == nn)]
    #else:
    #    df = df[(df["prev"] == triplet[0]) & (df["res"] == triplet[1]) & (df["next"] == triplet[2]) & (df["nn_pos"] == nn)]
    
    # Filter by phosphosite
    if phosphosite is not None:
        df = df[df["phos"] == phosphosite]

    return df[(df.prev == triplet[0]) & 
              (df.res == triplet[1]) & 
              (df.next == triplet[2]) & 
              (df.nn_res == nn)].protein_id.unique().tolist()

def get_sites_from_motif(
    df: pd.DataFrame,
    triplet: str, 
    nn: str, 
    ref_atom: str = None,
    phosphosite: bool = None,
) -> List[str]:
    """Return the list of sites (protein id, site) that contain a given motif."""
    
    # Filter by ref_atom
    if ref_atom is not None:
        if ref_atom == "CA": ref_atom = ["CA"]
        elif ref_atom == "oxygen": ref_atom = GAMMA_OXYGEN_CODES
        df = df[df["ref_atom"].isin(ref_atom)]
    
    # Turn dataframe cols into just residue (i.e. first character)
    for col in ["prev", "res", "next", "nn_res"]:
        df[col] = df[col].str[0] 

    # If "next" is unspecified: 
    #if len(triplet) == 2:
    #    df = df[(df["prev"] == triplet[0]) & (df["res"] == triplet[1]) & (df["nn_pos"] == nn)]
    #else:
    #    df = df[(df["prev"] == triplet[0]) & (df["res"] == triplet[1]) & (df["next"] == triplet[2]) & (df["nn_pos"] == nn)]
    
    # Filter by phosphosite
    if phosphosite is not None:
        df = df[df["phos"] == phosphosite]

    dff = df[(df.prev == triplet[0]) & 
              (df.res == triplet[1]) & 
              (df.next == triplet[2]) & 
              (df.nn_res == nn)]

    dff = dff[["protein_id", "res", "pos"]]
    
    # Return list of tuples (protein_id, res, pos)
    return list(dff.itertuples(index=False, name=None))

