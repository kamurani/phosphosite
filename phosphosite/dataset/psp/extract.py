# This code utilises the "PROTEINS API" supplied by UniProt. Within this is a section called proteomics-ptm.
# Given an accession number, a JSON string is returned containing all modifications for that
# protein.

# https://www.ebi.ac.uk/proteins/api/doc/#/

# This codes output a table with the following headings: 
# 'Accession', 'Entry', 'Mod Type', 'Mod Position', 'Mod Localization Probability'

from bioservices import UniProt
import requests, pprint
import pandas as pd
from bs4 import BeautifulSoup
import json

# Uniprot table extraction: Large-Scale data (almost all phosphorylation)
def get_large_scale_table(a):

    request_url = f"https://www.ebi.ac.uk/proteins/api/proteomics-ptm/{a}"
    r = requests.get(request_url, headers={ "Accept" : "application/json"})
    if not r.ok:
        return
    responseBody = r.json()

    # Get the relevant information from the JSON dictionary
    accession = responseBody['accession']
    entry = responseBody['entryName']
    modifications = responseBody['features']
    database = 'UniProt Large Scale'

    # Initialize lists to store the extracted information
    mod_type_list = []
    pep_start_list = []
    mod_pos_list = []
    prob_list = []
    
    # Loop over each modification and extract the relevant information
    for mod in modifications:
        mod_type_list.append(mod['ptms'][0]['name'])
        pep_start_list.append(mod['begin'])
        mod_pos_list.append(mod['ptms'][0]['position'])
        prob_list.append(mod['ptms'][0]['dbReferences'][0]['properties']['Localization probability'])

    # Finding true modification site
    # Position site by default is given by the position within the peptide
    # Mod site can be found by: peptide start site + position - 1
    for i in range(len(mod_pos_list)):
        mod_pos_list[i] = mod_pos_list[i] + int(pep_start_list[i]) - 1

    # Combine the extracted information into a Pandas DataFrame
    df = pd.DataFrame({
        'Accession': [accession] * len(modifications),
        'Entry': [entry] * len(modifications),
        'Mod Type': mod_type_list,
        'Mod Position': mod_pos_list,
        'Localization Probability': prob_list,
        'No of References': 'tbd1',
        'Database': [database] * len(modifications)
    })

    return(df)

# Finds the number of literature articles for a given modification
def count_references(references):
    count = 0
    for reference in references:
        try: 
            if reference['source']['name'] == 'PubMed':
                count += 1
        except KeyError:
            pass
    return count

# Uniprot table extraction: Non Large-Scale data
def get_non_large_scale_table(a):

    request_url = f"https://www.ebi.ac.uk/proteins/api/features/{a}"
    r = requests.get(request_url, headers={ "Accept" : "application/json"})
    if not r.ok:
        return
    responseBody = r.json()
    
    # Get the relevant information from the JSON dictionary
    modifications = responseBody['features']

    # Initialize lists to store the extracted information
    mod_type_list = []
    mod_pos_list = []
    reference_list = []

    # Loop over each modification and extract the relevant information
    mod_count = 0
    for mod in modifications:
        if mod['category'] == 'PTM' and mod['type'] == 'MOD_RES':
            mod_type_list.append(mod['description'])
            mod_pos_list.append(int(mod['begin']))
            reference_list.append(count_references(mod['evidences']))
            mod_count += 1
  
    # Combine the extracted information into a Pandas DataFrame
    df = pd.DataFrame({
        'Accession': responseBody['accession'],
        'Entry': responseBody['entryName'],
        'Mod Type': mod_type_list,
        'Mod Position': mod_pos_list,
        'Localization Probability': 'none',
        'No of References': reference_list,
        'Database': 'UniProt Non Large Scale'
    })

    return(df)

def get_manual_ptms(a):

    u = UniProt()
    res = u.retrieve(a, database='uniprot') # Database search
    
    accession = res['primaryAccession']
    entry = res['uniProtkbId']
    modifications = [x for x in res['comments'] if x.get('commentType') == 'PTM']
    database = 'UniProt Manual'
    mod_count = len(modifications)

    # Initialize lists to store the extracted information
    mod_type_list = []

    for dict_item in (res['comments']):
        if dict_item.get('commentType') == 'PTM':
            info, evidence = get_info_and_evidence(dict_item.get('texts'))
            mod_type_list.append(info)

    df = pd.DataFrame({
        'Accession': [accession] * mod_count,
        'Entry': [entry] * mod_count,
        'Mod Type': mod_type_list,
        'Mod Position': int(99999),
        'Localization Probability': 'none',
        'No of References': 'tbd3',
        'Database': [database] * mod_count
    })
    return df

# 'ptm_details' is a list of dicts. It always contains 'value' which describes the PTM.
# It can also have 'evidences' which cites relevant literature. Returns information and
# related evidence as a tuple.
def get_info_and_evidence(ptm_details):

    for dict_item in ptm_details:
        info = dict_item.get('value')
        if dict_item.get('evidences'):
            ev = get_evidence(dict_item.get('evidences'))
        else: 
            ev = 'NONE'
    return(info, ev)

def get_evidence(input):

    evidence = list()

    for dict_item in input:
        if dict_item.get('id') and dict_item.get('source'): 
            evidence.append(dict_item.get('source') + ' ' + dict_item.get('id'))
        elif dict_item.get('id') and not dict_item.get('source'):
            evidence.append(dict_item.get('id'))
        elif not dict_item.get('id') and dict_item.get('source'): 
            evidence.append(dict_item.get('source'))
        else:
            evidence.append('no source no id')
    return evidence


def get_psp_mods(protein, psp_id_dict):

    psp_id = psp_id_dict[protein][0]
    entry_name = psp_id_dict[protein][1]

    # MODIFICATION RETRIEVAL. For each protein, extract information on every modification found on that protein. 
    # Each modification is represented as a dictionary
    url = f"https://www.phosphosite.org/proteinAction?id={psp_id}&showAllSites=true" 
    response = requests.get(url) 
    html = response.content
    soup = BeautifulSoup(html, "html.parser") 

    ptms_param = soup.find('param', id = 'PTMsites') # PTMS are stored in the div class 'data container' under the param id 'PTMsites'. 
    if not ptms_param:
        print(f'No PTMS for this ribosomal protein {protein}')
        return
    ptms = ptms_param['value'] # This is a json string
    try:
        modifications = json.loads(ptms) # Converts json string to list of dicts
    except json.decoder.JSONDecodeError:
        return

    # DATA EXTRACTION AND TABLE CREATION. Relevant information about each modification is extracted and
    # outputted as a pandas data frame.
    mod_type_list = []
    mod_pos_list = []
    paper_num_list = []

    mod_count = 0
    for mod in modifications:
        mod_type_list.append(mod['MODIFICATION'])
        mod_pos_list.append(mod['POS'])
        paper_num_list.append(mod['REF'])
        mod_count += 1
  
    df = pd.DataFrame({
        'Accession': protein,
        'Entry': entry_name,
        'Mod Type': mod_type_list,
        'Mod Position': mod_pos_list,
        'Localization Probability': 'none',
        'No of References': paper_num_list,
        'Database': 'PSP' 
    })

    return df