from phosphosite import UNIPROT_SEQUENCE_PATH
from phosphosite.utils.io import load_seq_dict_from_file

sequence_dict = load_seq_dict_from_file(UNIPROT_SEQUENCE_PATH, key_format="uniprot_id")