"""Classes for interacting with datasets."""

from phosphosite.dataset.psp import PhosphoSequenceList, get_psp_regulatory_sites, phosphorylation, psp_filtered
from phosphosite.dataset.dbptm import dbptm_phosphorylation

__all__ = [
    "PhosphoSequenceList",
    "get_psp_regulatory_sites",
]

phosphorylation
psp_filtered

from phosphosite.dataset.uniprot import sequence_dict


