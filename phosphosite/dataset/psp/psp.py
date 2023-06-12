"""Classes for PhosphoSitePlus (PSP) data."""
# %%
# Phosphosite
# Author: Cam İmran <c.mcmenamie@unsw.edu.au>
# License: MIT
# Project Website: https://github.com/kamurani/phosphosite
# Code Repository: https://github.com/kamurani/phosphosite

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from typing import Union, List, Tuple, Dict
from pathlib import Path

import pandas as pd 
import numpy as np
import re

from phosphosite import PSP_DATASET_DIR



def load_dataset_file(
    filepath: Union[str, Path],
) -> pd.DataFrame:
    """Load the PhosphoSitePlus phosphorylation dataset."""
    if not Path(filepath).exists():
        raise ValueError(f"Filepath {filepath} does not exist.")
    
    
    compression = "gzip" if str(filepath).endswith(".gz") else None
    dataset = pd.read_csv(
        filepath,
        sep="\t",
        skiprows=3,
        compression=compression,
        low_memory=False,
    )
    return dataset

def get_phosphorylation_dataset(
    filepath: Union[str, Path] = PSP_DATASET_DIR / "Phosphorylation_site_dataset.gz",
    handle_isoforms: str = "remove", # "canonical", "keep"
) -> pd.DataFrame:
    """Load the PhosphoSitePlus phosphorylation dataset."""
    phosphorylation = load_dataset_file(filepath)

    if handle_isoforms not in ["keep", "remove", "canonical"]:
        raise ValueError(f"handle_isoforms must be one of ['keep', 'remove', 'canonical'], got {handle_isoforms}")
    if handle_isoforms == "canonical":
        raise NotImplementedError("canonical isoform handling not yet implemented")
    elif handle_isoforms == "keep":
        pass
    elif handle_isoforms == "remove":
        # Remove isoforms from phosphorylation dataset.
        # Turn {protein_id}-1 into {protein_id}
        old_total = len(phosphorylation[phosphorylation["ACC_ID"].str.contains("-")])
        # Count number of isoforms in ACC_ID column. 
        isoform1count = phosphorylation["ACC_ID"].str.split("-").str[1].value_counts().to_dict()[str(1)]

        # For all instances of {protein_id}-1 in ACC_ID column, 
        # replace with {protein_id}.
        phosphorylation["ACC_ID"] = phosphorylation["ACC_ID"].str.replace(
            re.compile(r"-1$"),
            "",
        )
        new_total = len(phosphorylation[phosphorylation["ACC_ID"].str.contains("-")])
        assert new_total == old_total - isoform1count

        # Now remove all rows containing dash. 
        phosphorylation = phosphorylation[~phosphorylation["ACC_ID"].str.contains("-")]

    return phosphorylation



"""Filter pd.Series of protein ids for isoforms."""
def filter_isoforms(
    method: str = "canonical", # "remove", "keep"

) -> pd.Series:
    """Filter pd.Series of protein ids for isoforms."""
    
    if method not in ["keep", "remove", "canonical"]:
        raise ValueError(f"method must be one of ['keep', 'remove', 'canonical'], got {method}")
    if method == "canonical":
        pass 

    # TODO
        


phosphorylation = get_phosphorylation_dataset()




"""
TODO: 
- this class should be more generalisable to other datasets. 
- have a way to instantiate this object with a different dataset.

"""

"""
TODO:
- create custom object for incorporating phosphosites (inherits seqrecord?)
- so can do seq.to_str() instead of ptm_seq.to_str(uniprot_id)
"""
class PhosphoSequenceList(object):
    """Stores a list of sequences from PhosphoSitePlus, with associated metadata.
    
    Sequence IDs are of the form GN:{gene_name}|{name}|{organism}|{uniprot_id} 
    Sequences include lowercase letters to indicate phosphorylation sites. 
    
    """

    def __init__(
        self, 
        fasta_path: Union[str, Path],
        verbose: bool = False,
        encoding: str = "ISO-8859-1", #'utf-8',
        file_format: str = "fasta",
        handle_isoforms: str = "keep", # "remove", "canonical"
        organism: str = None,

    ):
        """Initializes a PhosphoSequenceList object.

        This is a list of `SeqRecord`s loaded from PhosphoSitePlus, with associated metadata.

        Parameters
        ----------
        fasta_path : Union[str, Path]
            Path to file containing sequences from PhosphoSitePlus.
        verbose : bool, optional
            Whether to print verbose output, by default False
        encoding : str, optional
            Encoding of file, by default "ISO-8859-1"
        file_format : str, optional
            Format of file, by default "fasta"
        handle_isoforms : str, optional
            How to handle isoforms, by default "keep"
        organism : str, optional
            Organism to use, by default "human"
        
        Raises
        ------
        ValueError
            If file is not found.
        """
        self.fasta_path = Path(fasta_path)
        if not fasta_path.is_file():
            raise ValueError(f"File not found: {fasta_path}")
        self.verbose = verbose
        self.encoding = encoding
        self.file_format = file_format

        if handle_isoforms not in ["keep", "remove", "canonical"]:
            raise ValueError(f"handle_isoforms must be one of ['keep', 'remove', 'canonical'], got {handle_isoforms}")
        self.handle_isoforms = handle_isoforms
        self.organism = organism
        self._seqs = None

    @property
    def seqs(self):
        if self._seqs is None:
            self._seqs = self._get_seqs()
        return self._seqs
    
    def _get_seqs(self):
        """Load sequences from FASTA file."""
        seqs: List[SeqRecord] = self._load_seqs()
        n_original = len(seqs)
        new_seqs = []
        for seq in seqs:
            fields = seq.id.split("|")
            if len(fields) < 4:
                continue
            gene_name, name, organism, uniprot_id = fields

            # Check if sequence isoform 
            if "-" in uniprot_id:
                if self.handle_isoforms == "remove":
                    continue
                elif self.handle_isoforms == "canonical":
                    uniprot_id = uniprot_id.split("-")[0]

                    # TODO: map sequence as well somehow?
                elif self.handle_isoforms == "keep":
                    pass
                else:
                    raise ValueError(f"handle_isoforms must be one of ['keep', 'remove', 'canonical'], got {self.handle_isoforms}")
            # Filter by organism 
            if self.organism is not None and organism != self.organism:
                continue

            new_seqs.append(SeqRecord(
                seq=seq.seq,
                id=uniprot_id,
                annotations={
                    "gene_name": gene_name,
                    "organism": organism,
                },
                name=name,
                letter_annotations={
                    "phosphorylation_int": [
                        1 if c.islower() else 0
                        for c in str(seq.seq)
                    ],
                    "phosphorylation": "".join([
                        "P" if c.islower() else " "
                        for c in str(seq.seq)
                    ]),
                    # alternatively, put "p" or "-" 
                },
            ))
        if self.verbose: print(f"Using {len(new_seqs)} of {n_original} sequences from {self.fasta_path} (dropped {n_original - len(new_seqs)})")
        return new_seqs
    
    """
    Load from raw file.
    """
    def _load_seqs(
        self,
    ) -> List[SeqRecord]:
        """Load sequences from FASTA file.
        
        Returns
        -------
        List[SeqRecord]
            List of sequences.
        """
        with open(self.fasta_path, "r", encoding=self.encoding) as handle:
            seqs = list(SeqIO.parse(handle, self.file_format))
        if self.verbose: print(f"Loaded {len(seqs)} sequences from {self.fasta_path}")
        return seqs
    
    def _filter_seqs(
        self, 
        seqs: List[SeqRecord],

    ) -> List[SeqRecord]:
        """Filter a list of sequences."""

        pass

    """
    String representation of a given sequence and phosphorylation annotations.
    """
    def to_str(
        self,
        uniprot_id: str,
    ) -> str:
        """Returns the string representation of a given sequence and phosphorylation annotations.
        
        Parameters
        ----------
        uniprot_id : str
            UniProt ID of sequence.
            
        Returns
        -------
        str
            String representation of sequence and phosphorylation annotations.
        """
        seq = self.get_record(uniprot_id)

        # Get phosphorylation annotations
        phos = seq.letter_annotations["phosphorylation"]
        
        # Turn phos into string containing a "P" for every 1, and a " " for every 0
        phos_str = "".join(["P" if p else " " for p in phos])

        # Put phos underneath sequence
        return f"{seq.seq}\n{phos_str}"
    
    def __len__(self):
        return len(self.seqs)
    
    def __contains__(
        self, 
        i: str,
    ):
        for seq in self.seqs:
            if seq.id == i:
                return True
        return False
        
    
    def __getitem__(
        self, 
        i: Union[slice, str],
    ):
        if isinstance(i, slice):
            return self.seqs[i]
        elif isinstance(i, str):
            for seq in self.seqs:
                if seq.id == i:
                    return seq
            raise ValueError(f"Sequence not found: {i}")

        return self.seqs[i]
    
    def __iter__(self):
        return iter(self.seqs)
    
    def __repr__(self):
        return f"{self.__class__.__name__}({self.fasta_path})"
    
    def __str__(self):
        return f"{len(self)} sequences from {self.fasta_path}"

    """
    Get raw sequence as a string. 
    """
    def get_sequence(
        self, 
        i: Union[int, str],
    ) -> str: 
        """Returns string of Seq from SeqRecord."""
        return str(self.get_record(i).seq)

        
    def get_record(
        self, 
        i: Union[int, str],
    ) -> SeqRecord:
        """Returns the sequence with the given index or ID."""
        if isinstance(i, int):
            return self.seqs[i].seq
        elif isinstance(i, str):
            for seq in self.seqs:
                if seq.id == i:
                    return seq
            raise ValueError(f"Sequence not found: {i}")
    
    def get_id(self, i):
        return self.seqs[i].id
    
    def get_mod_ratio(
        self, 
        res_list: List[str],
    ) -> pd.DataFrame:
        res_list = [res.lower() for res in res_list]
        df = self.get_num_mods(res_list)
        for res in res_list:
            df[f"{res}_ratio"] = df[f"{res}_mod"] / df[f"{res}_all"]
        # Only keep ratio and id
        df = df[[f"{res}_ratio" for res in res_list] + ["id"]]
        # Express as log2 
        for res in res_list:
            df[f"{res}_ratio"] = -np.log2(df[f"{res}_ratio"])
        
        return df

    def get_num_mods(
        self, 
        res_list: List[str],
    ) -> pd.DataFrame:
        """Returns every sequence with number of modified residues and number of normal residues.

        Parameters
        ----------
            res_list (List[str]): List of modified residues to count.
        
        Returns
        -------
            pd.DataFrame: DataFrame with columns "id", "count", "residue", "type".
        
        """
        # Count number of lowercase and uppercase versions of residue, 
        # for each sequence.
        row_list = []
        for seq in self.seqs:
            seq_str = str(seq.seq)
            row = {
                "id": seq.id,
            }
            for res in res_list:
                row[f"{res.lower()}_mod"] = seq_str.count(res.lower())
                row[f"{res.lower()}_all"] = seq_str.count(res.lower()) + seq_str.count(res.upper())
            row_list.append(row)
        return pd.DataFrame(
            row_list,
            #columns=["id", "count", "residue", "type"],
        )
    
    def get_uniprot_ids(self):
        return [seq.id for seq in self.seqs]
    
"""
Get regulatory sites from PSP.
"""
def get_psp_regulatory_sites(
    path: Path,
    verbose: bool = False,
    handle_isoforms: str = "remove", # "keep", "remove", "canonical"
) -> pd.DataFrame:
    """Load the PhosphoSitePlus regulatory sites dataset.

    Parameters
    ----------
    path : Path
        Path to the PhosphoSitePlus regulatory sites dataset.
    verbose : bool
        Whether to print verbose output.
    handle_isoforms : str
        How to handle sequence isoforms.  One of "keep", "remove", "canonical".
        Default is "remove".

    Returns
    -------
    pd.DataFrame
        PhosphoSitePlus regulatory sites dataset.

    """
    isoform_options = ["keep", "remove", "canonical"]
    assert handle_isoforms in isoform_options, f"`handle_isoforms` must be one of {isoform_options}"
    
    df = pd.read_csv(
        path, 
        sep='\t',
        header=2,
    )
    df = df[df['ORGANISM'] == 'human']
    # Position, residue, and modification
    df[["RES", "POS", "MOD"]] = df["MOD_RSD"].str.extract('(?P<res>\w)(?P<pos>\d+)-(?P<mod>\w+)', expand=True)
    to_keep = ["PROTEIN", "ACC_ID", "ON_FUNCTION", "ON_PROCESS", "ON_PROT_INTERACT", "ON_OTHER_INTERACT",
               "RES", "POS", "MOD"]
    df = df[to_keep]
    # Deal with sequence isoforms
    if handle_isoforms == "remove":
        df = df[~df['ACC_ID'].str.contains('-')]    # Remove isoforms
    elif handle_isoforms == "canonical":
        df['ACC_ID'] = df['ACC_ID'].str.split('-').str[0] # Use canonical isoform only 
        # TODO: 
        # - other funcs will have to deal with this too in order to ensure that positions 
        # can be mapped onto sequences properly
    else:
        pass # Keep all isoforms as is

    return df

