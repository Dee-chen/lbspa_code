import os
import sys
import stereoAlign
import scanpy as sc
import numpy as np
from anndata import AnnData
import anndata as ad
sys.path.append(os.getcwd())

region = sys.argv[1]

lamprey = sc.read_h5ad(f'/data/work/22.fusemap/01.datas/{region}/lamrepy_brain_stereo_3d.h5ad')
mouse = sc.read_h5ad(f'/data/work/22.fusemap/01.datas/{region}/mouse_brain_stereo_3d.h5ad')
lamprey.var.index = [i.upper() for i in lamprey.var.index]
mouse.var.index = [i.upper() for i in mouse.var.index]


mouse.obs[['region', 'celltype']] = mouse.obs[['area_name', 'cell_subclass']].values
lamprey.obs[['region', 'celltype']] = lamprey.obs[['region', 'anno_whole_cluster']].values
mouse.obs['species'] = 'mouse'
lamprey.obs['species'] = 'lamprey'
adata = ad.concat([lamprey, mouse],)
adata.obs['species'] = adata.obs['species'].astype('category')



adata.obs.index = adata.obs.index.astype(str)
adata.obs_names_make_unique()
adata.layers["counts"] = adata.X

stereoAlign.pp.summarize_counts(adata)
stereoAlign.pp.norma_log(adata)
adata = stereoAlign.pp.scale_batch(adata, batch="species")


scvi_corrected = stereoAlign.alg.scvi_alignment(adata, batch_key="species", max_epochs = 200)
scvi_corrected.write(f'/data/work/22.fusemap/05.stereoalign/02.scvi/{region}/scvi_corrected.h5ad')


# python 05_crossSpecies_correlation.py cortex