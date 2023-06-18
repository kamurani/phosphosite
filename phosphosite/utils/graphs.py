"""Util functions for residue graphs."""

"""Sequence distance between two nodes."""
def get_seq_distance(
    node_id1: str, 
    node_id2: str,
    absolute: bool = True,
) -> int: 
    """Return difference in sequence position between two nodes.
    
    Parameters
    ----------
    node_id1 : str
        The first node ID.
    node_id2 : str
        The second node ID.
    abs : bool, optional
        Whether to return the absolute value, by default ``True``.

    Returns
    -------
    int
        The difference in sequence position between the two nodes.
    

    """
    pos1 = int(node_id1.split(":")[-1])
    pos2 = int(node_id2.split(":")[-1])
    diff = pos1 - pos2
    if absolute:
        return abs(diff)
    return diff