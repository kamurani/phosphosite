
import torch

from tqdm import tqdm 
from torch_geometric.data import Dataset, download_url
from torch_geometric.data import InMemoryDataset, download_url
from typing import List, Optional, Union, Tuple, Dict, Any, Callable
from pathlib import Path
from functools import partial
from graphein.protein.utils import download_alphafold_structure


# TODO: 
# raise issue in pytorch_geometric about error when creating dataset. (with download)
# Note: we had to resolve this issue manually by manually creating a `raw` directory 
# and moving files across that the `download` method was meant to create.

from phosphosite.graphs import get_pyg_graph

class PhosphoGraphDataset(InMemoryDataset):
    """
    Dataset for generating PyG graphs for phosphorylation sites
    on protein structures.

    Parameters
    ----------
    root : Path
        Path to directory containing raw and processed data.
    transform : Callable, optional
        Transform to apply to data, by default None
    pre_transform : Callable, optional
        Transform to apply to data before processing, by default None
    pre_filter : Callable, optional
        Filter to apply to data before processing, by default None

    """
    def __init__(
        self,
        root: Path,
        transform: Callable = None, 
        pre_transform: Callable = None, 
        pre_filter: Callable = None, 
        
        uniprot_ids: List[str] = None,

    ):
        self.uniprot_ids = uniprot_ids
        super().__init__(root, transform, pre_transform, pre_filter)
        self.data, self.slices = torch.load(self.processed_paths[0])
        #if uniprot_ids is None: 
        #    # Just use every file in the directory
        #    uniprot_ids = [f.stem for f in Path(self.raw_dir).glob('*.pdb')]  

        # TODO: this doesn't work because we need the 'raw dir' from the `super()`
        # call; but we can't call super() before we define the uniprot_ids. 
        #self.uniprot_ids = uniprot_ids
        

    @property
    def raw_file_names(self):
        # glob every .pdb file in the directory
        pdb_dir = self.raw_dir  
        return [f for f in pdb_dir.glob('*.pdb')] 

    @property
    def processed_file_names(self):
        return ['data.pt']
    

    #def download(self):
    #    for uniprot_id in self.uniprot_ids:
    #        download_alphafold_structure(uniprot_id, self.raw_dir)

    def process(self):
        
        func = partial(get_pyg_graph, pdb_dir=self.raw_dir)
        data_list = [func(uniprot_id) for uniprot_id in self.uniprot_ids]

        data_list = []
        successful = []
        pbar = tqdm(self.uniprot_ids)
        for i, uniprot_id in enumerate(pbar):
            pbar.set_description(f"Processing {uniprot_id}")
            pbar.update(1)
            try:
                g = func(uniprot_id)
            except Exception as e:
                print(f"Failed to generate graph for {uniprot_id}.")
                continue
            
            if g is not None:
                data_list.append(g)
                successful.append(uniprot_id)

    
        print(f"Successfully generated graphs for {len(successful)} out of {len(self.uniprot_ids)} proteins.")
        self.uniprot_ids = successful

        if self.pre_filter is not None:
            data_list = [data for data in data_list if self.pre_filter(data_list)]

        if self.pre_transform is not None:
            data_list = [self.pre_transform(data) for data in data_list]
        
        data, slices = self.collate(data_list)
        torch.save((data, slices), self.processed_paths[0])
        
from phosphosite import PHOSPHOSITE_PREDICT_DIR

#dataset = PhosphoGraphDataset(root=PHOSPHOSITE_PREDICT_DIR / "dataset")