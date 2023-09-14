
from typing import Tuple, List, Dict, Optional, Union


"""Split a string of the form "S119-p" into a tuple of the form `("S", "119", "p")`."""
def split_mod_rsd(mod_rsd: str) -> Tuple[str]:
    """Parse a MOD_RSD id into its components.

    Expects a string of the form "S119-p" and returns a tuple of the form `("S", 119)`.
    If there is no "-" character, assumes that the residue is unmodified and returns a 
    tuple of the form ("S", "119", None).
    
    
    """
    if "-" not in mod_rsd:
        # Split site into residue and position
        residue, position = mod_rsd[0], int(mod_rsd[1:])
        return residue, position
    # Remove modification after '-'
    site, mod = mod_rsd.split("-")
    # Split site into residue and position
    residue, position = site[0], int(site[1:])
    return residue, position, mod