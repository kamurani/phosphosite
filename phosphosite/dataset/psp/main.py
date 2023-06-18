import extract
import requests, pprint
import pandas as pd

###################
#### FUNCTIONS ####
###################

# Extracts modifications from UniProt and PSP and concatenates into a single dataframe
def get_ptms(a):

    df1 = extract.get_large_scale_table(a)
    df2 = extract.get_non_large_scale_table(a)
    df3 = extract.get_manual_ptms(a)
    psp_id_dict = get_psp_id_and_names()
    df4 = extract.get_psp_mods(a, psp_id_dict)
    df_all = pd.concat([df1, df2, df3, df4]).sort_values('Mod Position')
    return df_all

# Returns information about a protein's PSP ID and its UniProt entry name.
# This information can be used by using the UniProt accession.
def get_psp_id_and_names():

    df = pd.read_excel('../Human_Ribosome_List.xlsx')
    proteins = df.to_dict(orient='records')
    code_dict = dict()
    for protein in proteins:
        psp_id = str(protein['PSP ID']).strip()
        entry = protein['Entry'].strip()
        entry_name = protein['Entry Name'].strip()
        code_dict[entry] = [psp_id, entry_name]
    return code_dict

# All modification names are written into lower case. UniProt Manual entries are removed.
# Phosphotyrosine, phosphothreonine etc are re-written as 'phosphorylation'
# Any entries without a number in references is treated as 0 references
# Any entries without a localization probability is treated as 0 for the probability
def clean(df):
    df_copy = df.copy()
    df_copy['Mod Type'] = df_copy['Mod Type'].str.lower() 
    df_copy = df_copy[df_copy['Database'] != 'UniProt Manual'] 
    df_copy['Mod Type'] = df_copy['Mod Type'].replace(['phosphotyrosine', 'phosphoserine', 'phosphothreonine'], 'phosphorylation')
    df_copy['No of References'] = pd.to_numeric(df_copy['No of References'], errors='coerce').fillna(0) 
    df_copy['Localization Probability'] = pd.to_numeric(df_copy['Localization Probability'], errors='coerce').fillna(0)
    return df_copy 

# Some modifications will be found in multiple sources. In this case, all duplicates will be compressed into one row
# References are added, the highest localization probability value is selected, and a list of the databases is created 
def process_duplicates(df):
    groups = df.groupby(["Accession", "Entry", "Mod Type", "Mod Position"], sort=False)    
    result = groups.agg({"No of References": "sum", "Localization Probability": "max", "Database": lambda x: list(x)})
    filtered_df = result.reset_index()
    return filtered_df


###################
###### MAIN #######
###################

filepath = '../protein_list.txt'
with open(filepath, 'r') as f:
    protein_list = f.readlines()
    protein_list = [line.strip() for line in protein_list]
f.close

for protein in protein_list:
    df = get_ptms(protein)
    try:
        df_all = pd.concat([df_all, df])
    except NameError:
        df_all = df

df_all = clean(df_all)
filtered_df = process_duplicates(df_all)


pd.set_option('display.max_rows', None)   # Set option to display all rows
pd.set_option('display.max_columns', None)   # Set option to display all columns    
pd.set_option('display.width', None)
#print(df_all)

filtered_df.to_csv('modification_list_new.tsv', sep='\t', index=False)


