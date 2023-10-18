
import pandas as pd
from phosphosite import DBPTM_DATASET_DIR
phosphorylation_path = DBPTM_DATASET_DIR / "Phosphorylation"

def get_dbptm_phosphorylation():

    dbptm_phosphorylation = pd.read_csv(
        phosphorylation_path, 
        sep="\t",
        header=None,
    )
    dbptm_phosphorylation.columns = [
        "gene_name",
        "uniprot_id",
        "position",
        "ptm_type",
        "evidence",
        "seq_window",
    ]

    dbptm_phosphorylation["residue"] = dbptm_phosphorylation.apply(
        # Get middle residue of sequence window
        lambda x: x["seq_window"][int(len(x["seq_window"])/2)],
        axis=1,
    )

    # turn seq_window into str
    dbptm_phosphorylation["seq_window"] = dbptm_phosphorylation["seq_window"].astype(str)
    def get_middle_residue(seq_window):
        return seq_window[round(len(seq_window)/2)]

    dbptm_phosphorylation["residue"] = dbptm_phosphorylation["seq_window"].apply(get_middle_residue)    

    return dbptm_phosphorylation

dbptm_phosphorylation = get_dbptm_phosphorylation()