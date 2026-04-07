# =============================================================================
# 07_zebrafish_cerebellum.R — Zebrafish cerebellum spatial analysis
# =============================================================================
# Analyzes zebrafish cerebellum spatial transcriptomics data and compares
# cerebellar cell types with lamprey. Includes:
# - Spatial visualization of cerebellar layers (granular, Purkinje, molecular)
# - Cerebellar marker gene expression (pvalb7, ca8, aldoca, grid2)
# - Cross-species cerebellar region comparison with lamprey
#
# Input:  Zebrafish spatial data, lamprey-zebrafish integrated cerebellum object
# Output: Spatial plots, marker expression figures
# =============================================================================

source("00_utils.R")

# --- Load datasets -----------------------------------------------------------

# Zebrafish brain spatial transcriptomics
zebrafish.data <- readRDS(file.path(BASE_DIR,
  "20260324_zebrafishdata/adata.rds"))
zebrafish.data <- UpdateSeuratObject(zebrafish.data)
zebrafish.data <- NormalizeData(zebrafish.data)

# Lamprey-zebrafish integrated cerebellum object
cdy.cere <- readRDS(file.path(BASE_DIR,
  "20260324_zebrafishdata/01zebrafish_lamprey_spatial_geneIntersect_cb_analysis.rds"))

# Lamprey 3D brain atlas
data.lamprey.3D <- readRDS(LAMPREY_3D_NEW_PATH)


# --- Zebrafish cerebellar layer visualization ---------------------------------
# Highlight the four main cerebellar layers in spatial context

cerebellar_layers <- c("granular layer", "purkinje cell layer",
                       "molecular layer", "caudal lobe of cerebellum")
layer_colors <- c("#1d52a1", "#716db2", "#65c8cc", "#72c15a")

# Select a representative spatial section
section_names <- names(table(zebrafish.data$slice_code))

Idents(zebrafish.data) <- "region"

pdf(file.path(RESULT_DIR, "zebrafish_cerebellum.spatial.pdf"),
    width = 10, height = 8)

# Cerebellar layer spatial map
p <- Cluster_Highlight_Plot(
  seurat_object = subset(zebrafish.data,
                          slice_code == section_names[4]),
  group.by = "region", reduction = "spatial",
  cluster_name = cerebellar_layers,
  highlight_color = layer_colors)
print(p)

# Cerebellar marker gene expression in spatial context
cerebellar_markers <- c("pvalb7", "ca8", "aldoca", "grid2")
p2 <- FeaturePlot(
  subset(zebrafish.data, slice_code == section_names[4]),
  max.cutoff = 4, reduction = "spatial",
  features = cerebellar_markers, order = TRUE)
print(p2)

# GABAergic markers
gaba_markers <- c("slc32a1", "zic1", "gad1b", "gad2")
p3 <- FeaturePlot(
  subset(zebrafish.data, slice_code == section_names[4]),
  max.cutoff = 4, reduction = "spatial",
  features = gaba_markers, order = TRUE)
print(p3)

dev.off()


# --- Zebrafish cerebellum subset from integrated object -----------------------
# Map the integrated clustering back to zebrafish cerebellar cells

zebrafish.meta <- zebrafish.data@meta.data
zebrafish.meta.cere <- subset(zebrafish.meta,
  region %in% cerebellar_layers)

cdy.cere.zebra <- subset(cdy.cere, orig.ident == "sample")
cdy.cere.zebra <- subset(cdy.cere.zebra,
                           cells = rownames(zebrafish.meta.cere))

# Add anatomical region annotation from original zebrafish data
zebrafish.meta.cere <- zebrafish.meta.cere[
  rownames(cdy.cere.zebra@meta.data), ]
cdy.cere.zebra <- AddMetaData(cdy.cere.zebra,
  metadata = zebrafish.meta.cere[, 4], col.name = "region.new")
cdy.cere.zebra <- NormalizeData(cdy.cere.zebra)

pdf(file.path(RESULT_DIR, "zebrafish_cerebellum.integrated.pdf"),
    width = 10, height = 5)
# UMAP colored by anatomical region
p <- DimPlot(cdy.cere.zebra, group.by = "region.new", reduction = "umap",
             cols = c("#72c15a", "#1d52a1", "#65c8cc", "#716db2"))
print(p)
# Orthologous cerebellar markers in integrated space
p2 <- FeaturePlot(cdy.cere.zebra,
                   features = c("PVALB", "CA", "ALDOC", "GRID2"),
                   order = TRUE)
print(p2)
dev.off()


# --- Lamprey cerebellum and pons region spatial distribution ------------------

Idents(data.lamprey.3D) <- "region"

pdf(file.path(RESULT_DIR, "lamprey_cerebellum_pons.spatial.pdf"),
    width = 10, height = 8)
p <- Cluster_Highlight_Plot(
  seurat_object = subset(data.lamprey.3D, subset = slices == "28"),
  group.by = "region",
  cluster_name = c("cerebellum_area", "motor_nucleus_of_V",
                   "motor_nucleus_of_VII", "nucleus_of_X",
                   "areae_octavolateralis"),
  reduction = "align_spatial")
print(p)
dev.off()
