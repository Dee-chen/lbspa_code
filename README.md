# Lamprey 3D Single-Cell Transcriptomics Analysis Codebase

This repository contains the analysis code for the research article: "Lamprey 3D single-cell transcriptomics reveals ancestral and specialized features of the vertebrate brain". The codebase implements comprehensive spatial transcriptomics analyses of lamprey brain tissue, integrating multi-species comparisons to identify conserved and specialized vertebrate brain features.

## Repository Structure

### 1. `lamprey_spatial/`
Core spatial transcriptomics analysis pipeline for lamprey brain tissue:

| File/Directory | Description |
|----------------|-------------|
| `01_cell_segmentation/` | 40 Jupyter notebooks (`slice_1.ipynb` to `slice_40.ipynb`) containing cell segmentation workflows for individual brain slices |
| `02_spateo_align.ipynb` | Spatial alignment using Spateo framework for 3D reconstruction |
| `03_region_annotation.ipynb` | Brain region annotation and spatial clustering analysis |
| `04_graphst_clustering.py` | GraphST-based spatial clustering implementation |
| `05_crossSpecies_correlation.py` | Cross-species transcriptomic correlation analysis |
| `figs*.ipynb` | Figure-specific analysis notebooks for generating publication-quality visualizations (Figs 2, 8-16, 19, S2) |

### 2. `Brian_region&visualization_analysis/`
Cross-species comparative analysis and visualization workflows:

| Script | Description |
|--------|-------------|
| `00_utils.R` | Shared utility functions, package dependencies, and path configuration |
| `01_region_integration_analysis.R` | scVI-corrected data integration across 8 brain regions |
| `02_conserved_marker_analysis.R` | Identification of cross-species conserved marker genes |
| `03_spatial_visualization.R` | Spatial distribution plots and feature expression visualization |
| `04_pineal_cross_species.R` | Three-species pineal gland comparative analysis (lamprey, zebrafish, rat) |
| `05_non_neuron_analysis.R` | Non-neuronal cell type analysis (astroglia, ependymal cells) |
| `06_multi_species_neuron.R` | 8-vertebrate neuronal cell type comparison and gene specificity classification |
| `07_zebrafish_cerebellum.R` | Zebrafish cerebellum spatial transcriptomics analysis |

### 3. Top-Level Notebooks
| Notebook | Description |
|----------|-------------|
| `Single-cell_SpatialTranscriptomics_Integration_simplified.ipynb` | Harmony-based integration of single-cell and spatial transcriptomics data with batch correction |
| `Zebrafish _lamprey_Cerebellum_plot.ipynb` | Comparative cerebellum analysis between zebrafish and lamprey |

## Key Analyses Implemented

1. **3D Spatial Reconstruction**: Alignment and integration of 40 serial brain slices
2. **Cross-Species Comparisons**: Lamprey, mouse, zebrafish, rat, and 8-vertebrate panel analyses
3. **Cell Type Identification**: Comprehensive annotation of neuronal and non-neuronal cell populations
4. **Conserved Marker Discovery**: Evolutionarily conserved gene expression patterns across vertebrates
5. **Spatial Visualization**: Publication-quality 2D and 3D spatial transcriptomic plots

## Dependencies

- **R**: Seurat, dplyr, tidyverse, scCustomize, pheatmap, ggalluvial
- **Python**: scanpy, pandas, numpy, matplotlib, cv2, spateo, graphst
- **Single-Cell Tools**: scVI, Harmony, GraphST

## Usage

1. Clone the repository:
   ```bash
   git clone https://github.com/your-organization/lamprey-3d-brain.git