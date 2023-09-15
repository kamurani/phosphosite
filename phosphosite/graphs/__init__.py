"""CONSTANTS"""
NODE_DISTANCE_THRESHOLD = 6.0 # Å 
LONG_INTERACTION_THRESHOLD = 5 # 5 # How many sequence positions away can a node have its edges connected to it?


from phosphosite.graphs.subgraphs import get_motif_subgraph
from phosphosite.graphs.pyg import get_pyg_graph

get_motif_subgraph