"""Web scraping PhosphoSitePlus for phosphorylation data."""
# %%
# Phosphosite
# Author: Cam İmran <c.mcmenamie@unsw.edu.au>
# License: MIT
# Project Website: https://github.com/kamurani/phosphosite
# Code Repository: https://github.com/kamurani/phosphosite

# Modified from James Rosse <j.rosse@student.unsw.edu.au>
# Data Source: https://www.phosphosite.org 

import json
import requests
import pandas as pd 

from bs4 import BeautifulSoup

"""Extracts modifications from PhosphoSitePlus."""
def get_psp_mods(
    uniprot_id: str,
    base_url: str = "https://www.phosphosite.org/uniprotAccAction?id={uniprot_id}",
    verbose: bool = False,
) -> pd.DataFrame:
    """Loads the PhosphoSitePlus webpage and extracts the modifications 
    for a given protein.
    
    Parameters
    ----------
    uniprot_id : str
        The UniProt ID of the protein to extract modifications for.
    base_url : str
        The base URL to use for the PhosphoSitePlus website.
    
    Returns
    -------
    pd.DataFrame
        A pandas DataFrame containing the modifications for the given protein.

    """
    # MODIFICATION RETRIEVAL. 
    # For each protein, extract information on every modification found on that protein. 
    # Each modification is represented as a dictionary
    url = base_url.format(uniprot_id=uniprot_id)
    response = requests.get(url) 
    html = response.content
    soup = BeautifulSoup(html, "html.parser") 
    ptms_param = soup.find('param', id = 'PTMsites') # PTMS are stored in the div class 'data container' under the param id 'PTMsites'. 
    if not ptms_param:
        if verbose: print(f'No PTMS found for {uniprot_id}')
        return None
    ptms = ptms_param['value'] # This is a json string
    try:
        modifications = json.loads(ptms) # Converts json string to list of dicts
    except json.decoder.JSONDecodeError:
        raise ValueError(f"Could not decode JSON for {uniprot_id}")

    # DATA EXTRACTION AND TABLE CREATION. 
    # Relevant information about each modification is extracted and
    # formatted as a dataframe.
    df = pd.DataFrame.from_dict(modifications)
    return df
