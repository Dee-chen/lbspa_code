# Cross-species Brain Spatial Transcriptomics Analysis

Analysis code for the lamprey-mouse comparative brain spatial transcriptomics study. This directory contains the code contributed by one co-author, covering cross-species integration, conserved marker identification, and spatial visualization.

## Overview

This analysis compares cell types between **lamprey** (*Petromyzon marinus*) and **mouse** brain regions using spatial transcriptomics (Stereo-seq) integrated via scVI. Additional comparisons include zebrafish, rat, and a multi-species panel (8 vertebrates).

## File Structure

| Script | Description |
|--------|-------------|
| `00_utils.R` | Shared utility functions, package loading, and path configuration |
| `01_region_integration_analysis.R` | scVI-corrected data preprocessing and quality filtering for 8 brain regions |
| `02_conserved_marker_analysis.R` | Cross-species conserved marker identification at cluster and cell-class levels |
| `03_spatial_visualization.R` | Spatial distribution and feature expression plots for key brain regions |
| `04_pineal_cross_species.R` | Pineal gland three-species comparison (lamprey, zebrafish, rat) |
| `05_non_neuron_analysis.R` | Non-neuronal cell analysis (astroglia, ependymal cells) |
| `06_multi_species_neuron.R` | 8-species neuronal cell-type comparison and gene specificity classification |
| `07_zebrafish_cerebellum.R` | Zebrafish cerebellum spatial transcriptomics analysis |

## Brain Regions Covered

1. Olfactory bulb
2. Cortex (pallium)
3. Hippocampal formation
4. Striatum
5. Hypothalamus
6. Thalamus
7. Midbrain
8. Pons

## Dependencies

```r
# Core
Seurat (>= 4.0)
dplyr
tidyverse
patchwork
ggplot2

# Visualization
scCustomize
pheatmap
ggalluvial
RColorBrewer

# Statistics
matrixStats
```

## Usage

1. Edit `BASE_DIR` in `00_utils.R` to point to your data directory
2. Run scripts in numerical order (each script sources `00_utils.R`)
3. Ensure the `result/` output directory exists

## Data Availability

Input data (scVI-corrected Seurat objects, raw spatial stereo-seq data) are available from [TODO: add data repository link].
