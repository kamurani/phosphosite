"""Loading and annotation of structure files from AlphaFold and PDB."""

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
        filepath = self.get_filepath(uniprot_id)
        if not filepath.exists():
            raise ValueError(f"Filepath {filepath} does not exist.")
        return filepath