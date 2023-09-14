"""Loading and annotation of structure files from AlphaFold and PDB."""
import Bio.PDB
import gzip

from pathlib import Path
from typing import Union, List, Tuple, Dict, Optional

from phosphosite import AF_VERSION

# TODO: Retrieve uniprot_id from out_format f-string  
out_format = "AF-{uniprot_id}-F1-model_v{af_version}.{extension}"


class StructureLoader(object):
    af_version = AF_VERSION
    filename_template = "AF-{uniprot_id}-F1-model_v{af_version}.{extension}"

    def __init__(
        self,
        structure_dir: Path, 
        extension = "cif.gz",
        af_version = 3,
    ):
        self.structure_dir = structure_dir
        self.extension = extension
        self.af_version = af_version

    def get_filename(self, uniprot_id: str) -> str:
        return self.filename_template.format(
            uniprot_id=uniprot_id, 
            af_version=self.af_version,
            extension=self.extension,
        )
    
    def get_filepath(self, uniprot_id: str) -> Path:
        return self.structure_dir / self.get_filename(uniprot_id)
    
    def get_structure(self, uniprot_id: str) -> Path:
        """Return Path to structure file for given uniprot_id.
        
        Parameters
        ----------
        uniprot_id : str
            Uniprot ID of protein.
        
        Returns
        -------
        Path
            Path to structure file.

        Raises
        ------
        ValueError
            If structure file does not exist.
        """
        filepath = self.get_filepath(uniprot_id)
        if not filepath.exists():
            raise ValueError(f"Filepath {filepath} does not exist.")
        return filepath

    def protein_id_exists(self, uniprot_id: str) -> bool:
        """Check if structure file exists for given uniprot_id.

        Parameters
        ----------
        uniprot_id : str
            Uniprot ID of protein.
        
        Returns
        -------
        bool
            `True` if structure file exists, `False` otherwise.
        """
        return self.get_filepath(uniprot_id).exists()

    def get_existing_ids(
        self, 
        ids: Union[str, List[str]] = None,
    ) -> List[str]:
        """Returns intersection of given ids and existing ids."""
        if ids is not None:
            if isinstance(ids, str): ids = [ids]
            return [
                uniprot_id
                for uniprot_id in ids
                if self.protein_id_exists(uniprot_id)
            ]
        else:
            return [filepath.stem for filepath in self.structure_dir.glob(f"*.{self.extension}")]

    def parse_structure(
        self,
        uniprot_id: str,
    ) -> Bio.PDB.MMCIF2Dict.MMCIF2Dict:
        """Parse structure file for given uniprot_id.
        
        Parameters
        ----------
        uniprot_id : str
            Uniprot ID of protein.
        
        Returns
        -------
        Bio.PDB.MMCIF2Dict.MMCIF2Dict
            Structure file parsed as a dictionary.

        Raises
        ------
        ValueError
            If structure file has an unrecognized extension.
        
        """

        # get path to structure file.
        filepath: Path = self.get_structure(uniprot_id)
        if filepath.name.endswith("cif"):
            structure = Bio.PDB.MMCIF2Dict.MMCIF2Dict(filepath)
        elif filepath.name.endswith("cif.gz"):
            with gzip.open(filepath, "rt") as infile: 
                structure = Bio.PDB.MMCIF2Dict.MMCIF2Dict(infile)
        else: 
            raise ValueError(f"File '{filepath}' has an unrecognized extension.")

        return structure