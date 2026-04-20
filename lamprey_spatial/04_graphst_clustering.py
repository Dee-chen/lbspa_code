import sys
import torch
import scanpy as sc
from GraphST import GraphST


input_file = sys.argv[1]
output_file = sys.argv[2]

device = torch.device('cuda:3' if torch.cuda.is_available() else 'cpu')

adata = sc.read_h5ad(input_file)
adata.var_names_make_unique()

model = GraphST.GraphST(adata, device=device, random_seed=50)
adata = model.train()

adata.write(output_file)


# python 04_graphst_clustering.py rhombencephalon.h5ad rhombencephalon_graphst.h5ad