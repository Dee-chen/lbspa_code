# =============================================================================
# 00_utils.R — Shared utility functions and configuration
# =============================================================================
# This script provides common functions and path definitions used across all
# analysis scripts in this project. Source this file at the top of each script.
#
# Usage: source("00_utils.R")
# =============================================================================

# --- Required packages -------------------------------------------------------
library(Seurat)
library(dplyr)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(scCustomize)
library(pheatmap)
library(ggalluvial)
library(RColorBrewer)
library(matrixStats)

# --- Data path configuration -------------------------------------------------
# Adjust BASE_DIR to match your local directory structure.
# Expected layout:
#   BASE_DIR/
#   ├── scvi.raw/            # Raw scVI-corrected Seurat objects per region
#   │   ├── 01.olfactory_bulb/
#   │   ├── 02.cortex/
#   │   └── ...
#   ├── scvi_raw/            # Processed scVI objects (alternative naming)
#   ├── 01.datas/            # Raw spatial stereo-seq data per region
#   ├── lamprey-mouse/       # Cross-species reference data
#   ├── 17_pin/              # Pineal-specific datasets
#   └── result/              # Output directory for results

BASE_DIR    <- "."
SCVI_RAW    <- file.path(BASE_DIR, "scvi.raw")
SCVI_PROC   <- file.path(BASE_DIR, "scvi_raw")
SPATIAL_DIR <- file.path(BASE_DIR, "01.datas")
RESULT_DIR  <- file.path(BASE_DIR, "result")

# Path to the full lamprey 3D spatial object (used across multiple analyses)
LAMPREY_3D_PATH <- file.path(BASE_DIR, "lamprey.filter.fig1DE.rds")
LAMPREY_3D_NEW_PATH <- file.path(BASE_DIR, "lamprey.filter.fig1DE.newcluster.rds")

# Gene annotation file (BLAST NR results)
GENE_NR_PATH <- file.path(BASE_DIR, "blastp.NR.csv")

# --- Utility functions -------------------------------------------------------

#' Min-max scale a numeric vector to [0, 1]
#'
#' @param x Numeric vector
#' @return Scaled vector
min_max_scale <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

#' Create a custom spatial embedding from metadata coordinates
#'
#' Converts spatial coordinates stored in Seurat metadata into a
#' DimReducObject so that standard Seurat plotting functions can
#' render spatial layouts.
#'
#' @param seurat_obj Seurat object
#' @param coord_cols Character vector of length 2: column names in metadata
#'   for x and y coordinates
#' @param reduction_name Name for the new reduction (default: "custom_umap")
#' @return Modified Seurat object with the new reduction added
create_spatial_embedding <- function(seurat_obj, coord_cols = c("ax", "ay"),
                                     reduction_name = "custom_umap") {
  coords <- as.matrix(seurat_obj@meta.data[, coord_cols])
  colnames(coords) <- c("custom_1", "custom_2")
  seurat_obj[[reduction_name]] <- CreateDimReducObject(
    embeddings = coords,
    key = "custom_",
    assay = DefaultAssay(seurat_obj)
  )
  return(seurat_obj)
}

#' Run FindConservedMarkers across multiple clusters
#'
#' Iterates over a set of cluster identities and finds conserved markers
#' for each cluster between species.
#'
#' @param seurat_obj Seurat object with Idents set to cluster identity
#' @param clusters Vector of cluster IDs to iterate over
#' @param grouping_var Grouping variable for species (default: "species")
#' @param logfc_threshold Log-fold-change threshold (default: 0.1)
#' @param min_pct Minimum percentage threshold (default: 0.1)
#' @return Data frame with conserved markers for all clusters
run_conserved_marker_loop <- function(seurat_obj, clusters,
                                      grouping_var = "species",
                                      logfc_threshold = 0.1,
                                      min_pct = 0.1) {
  results <- list()
  for (cl in clusters) {
    tryCatch({
      markers <- FindConservedMarkers(
        seurat_obj,
        ident.1 = cl,
        logfc.threshold = logfc_threshold,
        min.pct = min_pct,
        grouping.var = grouping_var
      )
      markers$cluster <- cl
      results[[as.character(cl)]] <- markers
    }, error = function(e) {
      message(paste("Skipping cluster", cl, ":", e$message))
    })
  }
  do.call(rbind, results)
}

#' Standard Seurat preprocessing after scVI integration
#'
#' Runs NormalizeData, ScaleData, FindVariableFeatures, RunUMAP,
#' FindNeighbors, and FindClusters on scVI-aligned data.
#'
#' @param seurat_obj Seurat object with "aligned_scvi" reduction
#' @param dims Dimensions to use (default: 1:30)
#' @param resolutions Clustering resolutions (default: c(0.5, 1))
#' @return Processed Seurat object
preprocess_scvi_data <- function(seurat_obj, dims = 1:30,
                                  resolutions = c(0.5, 1)) {
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj)
  seurat_obj <- RunUMAP(seurat_obj, reduction = "aligned_scvi", dims = dims)
  seurat_obj <- FindNeighbors(seurat_obj, reduction = "aligned_scvi", dims = dims)
  for (res in resolutions) {
    seurat_obj <- FindClusters(seurat_obj, resolution = res)
  }
  return(seurat_obj)
}

#' Plot scaled heatmap of average gene expression
#'
#' Computes average expression per group, applies min-max scaling,
#' and generates a pheatmap.
#'
#' @param seurat_obj Seurat object
#' @param features Genes to plot
#' @param group_by Metadata column for grouping
#' @param color_palette Color for the heatmap gradient (paired with white)
#' @param title Plot title
#' @param cluster_rows Whether to cluster rows (default: FALSE)
plot_heatmap_scaled <- function(seurat_obj, features, group_by,
                                color_palette = "#C1640E",
                                title = "Marker Gene Expression (0-1 Scaled)",
                                cluster_rows = FALSE) {
  avg_exp <- AverageExpression(seurat_obj, features = features,
                               group.by = group_by, slot = "data")$RNA
  plot_matrix <- apply(avg_exp, 1, min_max_scale)
  pheatmap(plot_matrix,
           cluster_rows = cluster_rows,
           cluster_cols = FALSE,
           show_colnames = TRUE, angle_col = 45,
           scale = "none",
           color = colorRampPalette(c("white", color_palette))(100),
           border_color = "grey90",
           main = title)
}

#' Create a Sankey (alluvial) plot of cluster-to-celltype mapping
#'
#' Visualizes the correspondence between integrated clusters and
#' original cell type annotations using an alluvial diagram.
#'
#' @param meta_data Data frame (typically from seurat_obj@meta.data)
#' @param species_id Species identifier to filter (0 = lamprey, 1 = mouse)
#' @param celltype_indices Indices into the celltype factor levels to keep
#' @param cluster_col Column name for cluster assignment
#' @param celltype_col Column name for cell type annotation
#' @param freq_threshold Minimum frequency to include a connection
#' @param title Plot title
plot_sankey <- function(meta_data, species_id, celltype_indices,
                        cluster_col = "RNA_snn_res.0.5",
                        celltype_col = "celltype",
                        freq_threshold = 100,
                        title = "Cluster-CellType Sankey") {
  df <- subset(meta_data, species == species_id)
  df$cluster <- df[[cluster_col]]
  celltype_names <- names(table(df[[celltype_col]]))[celltype_indices]
  df <- subset(df, get(celltype_col) %in% celltype_names)

  plot_data <- df %>% count(cluster, !!sym(celltype_col), name = "Freq")

  filtered_data <- plot_data %>% filter(Freq >= freq_threshold)

  # Order celltypes by their dominant cluster
  cluster_levels <- unique(filtered_data$cluster)
  celltype_order <- filtered_data %>%
    group_by(!!sym(celltype_col)) %>%
    summarise(total_freq = sum(Freq),
              main_cluster = cluster[which.max(Freq)]) %>%
    arrange(match(main_cluster, cluster_levels), total_freq) %>%
    pull(!!sym(celltype_col))

  filtered_data <- filtered_data %>%
    mutate(
      cluster = factor(cluster, levels = cluster_levels),
      !!celltype_col := factor(get(celltype_col), levels = celltype_order)
    )

  ggplot(filtered_data,
         aes(y = Freq, axis1 = cluster, axis2 = !!sym(celltype_col))) +
    geom_alluvium(aes(fill = !!sym(celltype_col)),
                  width = 1/12, alpha = 0.7) +
    geom_stratum(width = 1/12, fill = "lightgray", color = "white") +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
    scale_x_discrete(limits = c("Cluster", "CellType"),
                     expand = c(0.05, 0.05)) +
    labs(title = title) +
    theme_minimal()
}
