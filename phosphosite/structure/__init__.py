
import functools
from phosphosite.structure.loader import StructureLoader
from phosphosite import AF_HUMAN_CIF as structure_dir

"""Default structure loader. 

Uses the user-specified structure directory (AF_HUMAN_CIF) from 
phosphosite/notebooks/config.yml 
"""
loader = StructureLoader(structure_dir=structure_dir) 

# Functools partial to make it easier to use the loader.
get_structure = functools.partial(loader.get_structure)

