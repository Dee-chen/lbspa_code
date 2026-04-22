# =============================================================================
# 03_spatial_visualization.R — Spatial transcriptomics visualization
# =============================================================================
# Generates key spatial visualizations for each brain region, including:
# - Spatial cell class distribution on tissue sections
# - Spatial feature plots for conserved marker genes
# - Sankey (alluvial) plots for cluster-celltype correspondence
#
# Input:  Filtered Seurat objects + raw spatial stereo-seq data
# Output: PDF figures for manuscript
# =============================================================================

source("00_utils.R")

# =============================================================================
# SECTION 1: Olfactory bulb — spatial cell class + feature plots
# =============================================================================

data <- readRDS(file.path(RESULT_DIR, "01.olfactory_bulb.rds"))
data <- NormalizeData(data)

# --- Subset and assign classes (same mapping as 02_conserved_marker_analysis) ---
ct_names <- names(table(data$celltype))
data.mouse   <- subset(data, subset = species == 1 &
                         celltype %in% ct_names[c(1, 6, 7, 34:37, 67:69)])
data.lamprey <- subset(data, subset = species == 0 &
                         celltype %in% ct_names[16:33])

data.mouse$class   <- data.mouse$celltype
levels(data.mouse$class)   <- c("Astro", "Endo", "Epend", "Exc", "Inh",
                                 "OPC", "Oligo", "Exc", "Inh", "Endo")
data.lamprey$class <- data.lamprey$celltype
levels(data.lamprey$class) <- c("Glia-like", "Glia-like", "Glia-like",
                                 "Inh-like", "Glia-like", "Glia-like",
                                 "Glia-like", "Inh-like", "Exc-like",
                                 "Inh-like", "Inh-like", "Exc-like",
                                 "Inh-like", "Exc-like", "Glia-like",
                                 "Glia-like", "Inh-like", "Glia-like")

# Create spatial embeddings from stereo-seq coordinates
data.mouse   <- create_spatial_embedding(data.mouse, c("ax", "ay"))
data.lamprey <- create_spatial_embedding(data.lamprey, c("ax", "ay"))

# --- Mouse spatial plots (representative sections T253, T257) ---
mouse_cols <- c("#F7DB9A", "#909090", "#D39D6E", "#FA6868", "#489EE9",
                "#5DD39E", "#24905F")
ob_genes <- c("SPARC", "ATP1B2", "ADGRB1", "GRIA3", "PRKCE",
              "GAD1", "MEIS2", "SLC32A1", "KCND2")

pdf(file.path(RESULT_DIR, "01.olfactory_bulb.spatial.pdf"),
    width = 10, height = 8)
# Cell class distribution on tissue
for (sl in c("T253", "T257")) {
  p <- DimPlot(subset(data.mouse, subset = slices == sl),
               group.by = "class", reduction = "custom_umap",
               cols = mouse_cols) + ggtitle(paste("Mouse OB -", sl))
  print(p)
}
# Conserved gene expression in spatial context
for (sl in c("T253", "T257")) {
  p <- FeaturePlot_scCustom(
    seurat_object = subset(data.mouse, subset = slices == sl),
    reduction = "custom_umap", max.cutoff = 4, features = ob_genes)
  print(p)
}
dev.off()

# --- Lamprey spatial plots (representative sections 3, 6) ---
lamprey_cols <- c("#909090", "#5A9CB5", "#FA6868")

pdf(file.path(RESULT_DIR, "01.olfactory_bulb.lamprey_spatial.pdf"),
    width = 10, height = 8)
for (sl in c("3", "6")) {
  p <- DimPlot(subset(data.lamprey, subset = az == sl),
               group.by = "class", reduction = "custom_umap",
               cols = lamprey_cols) + ggtitle(paste("Lamprey OB - section", sl))
  print(p)
  p2 <- FeaturePlot_scCustom(
    seurat_object = subset(data.lamprey, subset = az == sl),
    reduction = "custom_umap", max.cutoff = 4, features = ob_genes)
  print(p2)
}
dev.off()

# --- Sankey plots for olfactory bulb ---
pdf(file.path(RESULT_DIR, "01.olfactory_bulb.sankey.pdf"),
    width = 10, height = 8)
plot_sankey(data@meta.data, species_id = 0, celltype_indices = 16:33,
            freq_threshold = 282, title = "Olfactory Bulb - Lamprey")
plot_sankey(data@meta.data, species_id = 1,
            celltype_indices = c(1, 6, 7, 34:37, 67:69),
            freq_threshold = 577, title = "Olfactory Bulb - Mouse")
dev.off()

rm(data, data.mouse, data.lamprey)


# =============================================================================
# SECTION 2: Cortex — spatial cell class + layer marker visualization
# =============================================================================

data <- readRDS(file.path(SCVI_PROC, "02.cortex.rds"))
data <- NormalizeData(data)

ct_names <- names(table(data$celltype))
data.mouse   <- subset(data, subset = species == 1 &
                         celltype %in% ct_names[c(1:11, 31:33)])
data.lamprey <- subset(data, subset = species == 0 &
                         celltype %in% ct_names[12:30])

# Class mapping
data.mouse$class   <- data.mouse$celltype
levels(data.mouse$class) <- c("Astro", "Endo", "Gran-DG", "Endo", "Epend",
                               "Epend", "Mic", "Exc", "Inh", "OPC", "Oligo",
                               "Exc", "Inh", "Endo")
data.lamprey$class <- data.lamprey$celltype
levels(data.lamprey$class) <- c("Astro-like", "Exc-like2", "Exc-like2",
                                 "Astro-like", "Exc-like1", "Exc-like3",
                                 "Exc-like3", "Exc-like3", "Epend-like",
                                 "Astro-like", "Exc-like2", "Exc-like3",
                                 "Exc-like1", "Exc-like2", "Inh-like",
                                 "Inh-like", "Inh-like", "Astro-like",
                                 "Epend-like")

# Spatial embeddings
data.mouse$tx <- as.numeric(as.character(data.mouse$x))
data.mouse$ty <- as.numeric(as.character(data.mouse$y))
data.mouse <- create_spatial_embedding(data.mouse, c("ty", "tx"))

data.lamprey$x <- as.numeric(data.lamprey$x)
data.lamprey$y <- as.numeric(data.lamprey$y)
data.lamprey <- create_spatial_embedding(data.lamprey, c("x", "y"))

# Cortical layer marker genes
layer_genes <- c("RASGRF2", "VIP", "CARTPT", "KCNIP2", "MYL4", "RAB3C",
                 "RPRM", "CRYAB", "IGSF21", "TBR1", "CUX1")
cortex_genes <- c("SLC1A2", "PHGDH", "KCNIP4", "VSNL1", "DPP10",
                  "GAD1", "KCNC1", "CPLX1")

# Mouse cortex layers
data.mouse$layer <- str_extract(data.mouse$region, "\\d+$")

pdf(file.path(RESULT_DIR, "02.cortex.spatial.pdf"), width = 12, height = 8)
# Mouse sections T307, T326
for (sl in c("T307", "T326")) {
  p <- DimPlot(subset(data.mouse, subset = slices == sl),
               group.by = "class", reduction = "custom_umap") +
    ggtitle(paste("Mouse Cortex -", sl))
  print(p)
  p2 <- FeaturePlot_scCustom(
    seurat_object = subset(data.mouse, subset = slices == sl),
    reduction = "custom_umap", max.cutoff = 5, order = TRUE,
    features = cortex_genes)
  print(p2)
  # Layer markers
  p3 <- FeaturePlot_scCustom(
    seurat_object = subset(data.mouse, subset = slices == sl),
    reduction = "custom_umap", max.cutoff = 5, order = TRUE,
    features = layer_genes)
  print(p3)
}
# Lamprey sections 6, 35
lamprey_ctx_cols <- c("#909090", "#FA6868", "#F99569", "#E5D613",
                      "#5E5B25", "#489EE9")
for (sl in c("6", "35")) {
  p <- DimPlot(subset(data.lamprey, subset = slices == sl),
               group.by = "class", reduction = "custom_umap",
               cols = lamprey_ctx_cols) +
    ggtitle(paste("Lamprey Cortex - section", sl))
  print(p)
  p2 <- FeaturePlot_scCustom(
    seurat_object = subset(data.lamprey, subset = slices == sl),
    reduction = "custom_umap", max.cutoff = 5, features = cortex_genes)
  print(p2)
}
dev.off()

rm(data, data.mouse, data.lamprey)


# =============================================================================
# SECTION 3: Hypothalamus — spatial visualization
# =============================================================================

data <- readRDS(file.path(SCVI_PROC, "05.Hypothalamus.rds"))
data <- NormalizeData(data)

ct_names <- names(table(data$celltype))
data.mouse   <- subset(data, subset = species == 1 &
                         celltype %in% ct_names[c(1:6, 27:30, 33)])
data.lamprey <- subset(data, subset = species == 0 &
                         celltype %in% ct_names[7:26])

# Class mapping
data.mouse$class   <- data.mouse$celltype
levels(data.mouse$class) <- c("Astro", "Mono", "Chor", "DiMes", "Endo",
                               "Epend", "Epend", "Micro", "OPC", "Oligo", "Endo")
data.lamprey$class <- data.lamprey$celltype
levels(data.lamprey$class) <- c("Glia-like", "Glia-like", "Glia-like",
                                 "DiMes-like", "DiMes-like", "Mono-like",
                                 "DiMes-like", "Mono-like", "DiMes-like",
                                 "DiMes-like", "HypNeurons", "Glia-like",
                                 "HypNeurons", "Mono-like", "DiMes-like",
                                 "Epend-like", "Glia-like", "Glia-like",
                                 "Epend-like", "Epend-like")

# Spatial embeddings using raw 3D stereo-seq coordinates
data.mouse.raw <- readRDS(file.path(SPATIAL_DIR,
                                     "05.Hypothalamus/mouse_brain_stereo_3d.rds"))
data.mouse.raw <- subset(data.mouse.raw, cells = colnames(data.mouse))
data.mouse <- create_spatial_embedding(data.mouse,
  coord_cols = colnames(data.mouse.raw@meta.data[, c("ax", "ay")]))

data.lamprey.raw <- readRDS(file.path(SPATIAL_DIR,
                                       "05.Hypothalamus/lamrepy_brain_stereo_3d.rds"))
data.lamprey.raw <- subset(data.lamprey.raw, cells = colnames(data.lamprey))
data.lamprey <- create_spatial_embedding(data.lamprey, c("x", "y"))

hyp_genes <- c("CHGA", "CPLX1", "PCDH9", "SLC17A6", "GAD1", "CNP")

pdf(file.path(RESULT_DIR, "05.hypothalamus.spatial.pdf"), width = 10, height = 8)
# Mouse section T317
p <- DimPlot(subset(data.mouse, subset = slices == "T317"),
             group.by = "class", reduction = "custom_umap",
             cols = c("#909090", "#FA6868", "#5D5D93", "#F99569",
                      "#4C481D", "#4C1D2C", "#2D2D2D", "#5DD39E", "#24905F"))
print(p)
p2 <- FeaturePlot_scCustom(
  seurat_object = subset(data.mouse, subset = slices == "T317"),
  reduction = "custom_umap", max.cutoff = 4, features = hyp_genes)
print(p2)

# Lamprey section 25
p <- DimPlot(subset(data.lamprey, subset = slices == "25"),
             group.by = "class", reduction = "custom_umap",
             cols = c("#909090", "#5A9CB5", "#FA6868", "#F99569", "grey"))
print(p)
p2 <- FeaturePlot_scCustom(
  seurat_object = subset(data.lamprey, subset = slices == 25),
  reduction = "custom_umap", max.cutoff = 4, order = TRUE, alpha_exp = 0.5,
  features = hyp_genes)
print(p2)
dev.off()

rm(data, data.mouse, data.lamprey, data.mouse.raw, data.lamprey.raw)


# =============================================================================
# SECTION 4: Midbrain — spatial cell class + monoamine gene expression
# =============================================================================

data <- readRDS(file.path(SCVI_PROC, "07.Midbrain.rds"))
data <- NormalizeData(data)

ct_names <- names(table(data$celltype))
data.mouse   <- subset(data, subset = species == 1 &
                         celltype %in% ct_names[c(1:3, 5:9, 32:34, 36)])
data.lamprey <- subset(data, subset = species == 0 &
                         celltype %in% ct_names[10:30])

# Class mapping
data.mouse$class   <- data.mouse$celltype
levels(data.mouse$class) <- c("Astro", "Mono", "Chor", "DiMes", "Endo",
                               "Epend", "Epend", "Mic", "OPC", "Oligo",
                               "Rh", "Endo")
data.lamprey$class <- data.lamprey$celltype
levels(data.lamprey$class) <- c("Glia-like1", "MidNeuron", "MidNeuron",
                                 "MidNeuron", "DiMes-like1", "Mono-like",
                                 "Rh-like", "MidNeuron", "Rh-like",
                                 "DiMes-like1", "DiMes-like1", "Mono-like",
                                 "Rh-like", "Mono-like", "Mono-like",
                                 "DiMes-like1", "MidNeuron", "MidNeuron",
                                 "Glia-like2", "Glia-like1", "Glia-like2")

# Load raw spatial data for coordinate-based plotting
data.mouse.raw <- readRDS(file.path(SPATIAL_DIR,
                                     "07.Midbrain/mouse_brain_stereo_3d.rds"))
data.mouse.raw <- subset(data.mouse.raw,
  subset = cell_subclass %in% names(table(data.mouse.raw$cell_subclass))[
    c(1:3, 5:9, 11:13, 15)])
data.mouse.raw$class <- data.mouse.raw$cell_subclass
levels(data.mouse.raw$class) <- c("Astro", "Mono", "Chor", "DiMes", "Endo",
                                   "Epend", "Epend", "Mic", "OPC", "Oligo",
                                   "Rh", "Endo")
data.mouse.raw <- NormalizeData(data.mouse.raw)
data.mouse.raw <- create_spatial_embedding(data.mouse.raw, c("ax", "ay"))

data.lamprey <- create_spatial_embedding(data.lamprey, c("ax", "ay"))

# Conserved structural genes + monoamine pathway genes
midbrain_struct_genes <- c("ALDOC", "RTN1", "NEFL", "GAD1", "SLC17A6",
                           "GRIA2", "CNTNAP2", "SPARC")
monoamine_genes <- c("TH", "SLC18A2", "TPH2", "SLC6A2", "SNCG", "ADRB1")

pdf(file.path(RESULT_DIR, "07.midbrain.spatial.pdf"), width = 12, height = 10)

# Mouse midbrain sections T334, T341
for (sl in c("T334", "T341")) {
  # Cell class highlight plots
  p <- Cluster_Highlight_Plot(
    seurat_object = subset(data.mouse.raw, subset = slices == sl),
    group.by = "class",
    cluster_name = c("Mono", "DiMes", "Rh"),
    highlight_color = c("#FA6868", "#F99569", "#A35BB5"))
  print(p)
  # Structural gene expression
  p2 <- FeaturePlot_scCustom(
    seurat_object = subset(data.mouse.raw, subset = slices == sl),
    reduction = "custom_umap", max.cutoff = 5,
    features = str_to_title(midbrain_struct_genes))
  print(p2)
  # Monoamine biosynthesis genes
  p3 <- FeaturePlot_scCustom(
    seurat_object = subset(data.mouse.raw, subset = slices == sl),
    reduction = "custom_umap", max.cutoff = 5,
    features = str_to_title(monoamine_genes))
  print(p3)
}

# Lamprey midbrain sections 14, 26
for (sl in c("14", "26")) {
  p <- Cluster_Highlight_Plot(
    seurat_object = subset(data.lamprey, subset = slices == sl),
    group.by = "class",
    cluster_name = c("Mono-like", "DiMes-like1", "Rh-like", "MidNeuron"),
    highlight_color = c("#FA6868", "#F99569", "#A35BB5", "#489EE9"))
  print(p)
  p2 <- FeaturePlot_scCustom(
    seurat_object = subset(data.lamprey, subset = slices == sl),
    reduction = "custom_umap", max.cutoff = 6,
    features = midbrain_struct_genes)
  print(p2)
  p3 <- FeaturePlot_scCustom(
    seurat_object = subset(data.lamprey, subset = slices == sl),
    reduction = "custom_umap", max.cutoff = 6,
    features = monoamine_genes)
  print(p3)
}
dev.off()

# --- Midbrain Sankey plots ---
pdf(file.path(RESULT_DIR, "07.midbrain.sankey.pdf"), width = 10, height = 8)
plot_sankey(data@meta.data, species_id = 0, celltype_indices = 7:27,
            freq_threshold = 323, title = "Midbrain - Lamprey")
plot_sankey(data@meta.data, species_id = 1,
            celltype_indices = c(1:6, 28:31, 35),
            freq_threshold = 285, title = "Midbrain - Mouse")
dev.off()

rm(data, data.mouse, data.lamprey, data.mouse.raw)


# =============================================================================
# SECTION 5: Pons �� spatial cell class + locus coeruleus analysis
# =============================================================================

data <- readRDS(file.path(SCVI_PROC, "08.Pons.rds"))
data <- NormalizeData(data)

ct_names <- names(table(data$celltype))
data.mouse   <- subset(data, subset = species == 1 &
                         celltype %in% ct_names[c(1, 17:28, 50:53)])
data.lamprey <- subset(data, subset = species == 0 &
                         celltype %in% ct_names[2:16])

# Class mapping
data.mouse$class   <- data.mouse$celltype
levels(data.mouse$class) <- c("Astro", "Neuroblast", "Mono", "Chor", "DiMes",
                               "Endo", "Epend", "Epend", "Mic", "OPC", "Oligo",
                               "Rh", "Exc/Inh", "Exc/Inh", "Endo")
data.lamprey$class <- data.lamprey$celltype
levels(data.lamprey$class) <- c("Astro-like", "Epend-like", "Epend-like",
                                 "Epend-like", "Exc/Inh-like", "Neuroblast",
                                 "Rh-like", "Neuroblast", "Neuroblast",
                                 "Neuroblast", "Mono-like", "Rh-like",
                                 "Mono-like", "Mono-like", "Exc/Inh-like",
                                 "Rh-like", "Epend-like", "Astro-like",
                                 "Exc/Inh-like", "Mono-like")

data.mouse   <- create_spatial_embedding(data.mouse, c("ax", "ay"))
data.lamprey <- create_spatial_embedding(data.lamprey, c("ax", "ay"))

# Locus coeruleus (LC) marker comparison: mouse LP-Psat-LC vs lamprey Mono-like
data.mouse.left <- subset(data.mouse,
  subset = region %in% names(table(data.mouse$region))[1:94])
Idents(data.mouse.left) <- "region"
LC_markers <- FindMarkers(data.mouse.left, ident.1 = "LP-Psat-LC",
                          logfc.threshold = 0.5, min.pct = 0.1)
Idents(data.lamprey) <- "class"
Mono_markers <- FindMarkers(data.lamprey, ident.1 = "Mono-like",
                            logfc.threshold = 0.25, min.pct = 0.1)
LC_conserved <- intersect(rownames(LC_markers), rownames(Mono_markers))
message("LC-Mono conserved markers: ", length(LC_conserved))

pdf(file.path(RESULT_DIR, "08.pons.spatial.pdf"), width = 12, height = 10)
# Mouse pons sections T349, T352
mouse_pons_cols <- c("#909090", "#FA6868", "#5D5D93", "#F99569", "#4C481D",
                     "#4C1D2C", "#2D2D2D", "#5DD39E", "#24905F", "#A35BB5")
for (sl in c("T349", "T352")) {
  p <- DimPlot(subset(data.mouse.left, subset = slices == sl),
               group.by = "class", reduction = "custom_umap",
               cols = mouse_pons_cols)
  print(p)
  # Highlight LC region
  p2 <- Cluster_Highlight_Plot(
    seurat_object = subset(data.mouse.left, subset = slices == sl),
    cluster_name = c("LP-Psat-LC"), highlight_color = c("#A35BB5"),
    reduction = "custom_umap")
  print(p2)
}
# Lamprey pons section 20
p <- Cluster_Highlight_Plot(
  seurat_object = subset(data.lamprey, subset = slices == "20"),
  group.by = "class",
  cluster_name = c("Mono-like"), highlight_color = c("#A35BB5"),
  reduction = "custom_umap")
print(p)
# LC conserved gene expression in lamprey
if (length(LC_conserved) > 0) {
  p2 <- FeaturePlot_scCustom(
    seurat_object = subset(data.lamprey, subset = slices == "20"),
    reduction = "custom_umap", max.cutoff = 5, features = LC_conserved)
  print(p2)
}
dev.off()

rm(data, data.mouse, data.lamprey, data.mouse.left)


# =============================================================================
# SECTION 6: Thalamus — spatial visualization
# =============================================================================

data <- readRDS(file.path(SCVI_PROC, "06.Thalamus.rds"))
data <- NormalizeData(data)

ct_names <- names(table(data$celltype))
data.mouse   <- subset(data, subset = species == 1 &
                         celltype %in% ct_names[c(1:8, 10, 11, 32)])
data.lamprey <- subset(data, subset = species == 0 &
                         celltype %in% ct_names[15:31])

data.mouse$class   <- data.mouse$celltype
levels(data.mouse$class) <- c("Astro", "Mono", "Chor", "DiMes", "Endo",
                               "Epend", "Epend", "Micro", "OPC", "Oligo", "Endo")
data.lamprey$class <- data.lamprey$celltype
levels(data.lamprey$class) <- c("Glia-like", "Glia-like", "DiMes-like",
                                 "Mono-like", "Mono-like", "DiMes-like",
                                 "Endo-like", "Endo-like", "Endo-like",
                                 "Mono-like", "DiMes-like", "DiMes-like",
                                 "Mono-like", "Mono-like", "DiMes-like",
                                 "Mono-like", "Glia-like")

# Spatial embeddings from raw data
data.mouse.raw <- readRDS(file.path(SPATIAL_DIR,
                                     "06.Thalamus/mouse_brain_stereo_3d.rds"))
data.mouse.raw <- subset(data.mouse.raw, cells = colnames(data.mouse))
data.mouse <- create_spatial_embedding(data.mouse,
  coord_cols = c("ax", "ay"))

data.lamprey.raw <- readRDS(file.path(SPATIAL_DIR,
                                       "06.Thalamus/lamrepy_brain_stereo_3d.rds"))
data.lamprey.raw <- subset(data.lamprey.raw, cells = colnames(data.lamprey))
data.lamprey <- create_spatial_embedding(data.lamprey, c("x", "y"))

thal_genes <- c("SLC17A6", "GAD1", "CBLN1", "ZIC1", "NTNG1")

pdf(file.path(RESULT_DIR, "06.thalamus.spatial.pdf"), width = 10, height = 8)
# Mouse section T307
p <- DimPlot(subset(data.mouse, subset = slices == "T307"),
             group.by = "class", reduction = "custom_umap")
print(p)
p2 <- FeaturePlot_scCustom(
  seurat_object = subset(data.mouse, subset = slices == "T307"),
  reduction = "custom_umap", max.cutoff = 4, features = thal_genes)
print(p2)
# Lamprey section 25
p <- DimPlot(subset(data.lamprey, subset = slices == "25"),
             group.by = "class", reduction = "custom_umap",
             cols = c("#909090", "#5A9CB5", "#FA6868", "grey"))
print(p)
p2 <- FeaturePlot_scCustom(
  seurat_object = subset(data.lamprey, subset = slices == 25),
  reduction = "custom_umap", max.cutoff = 4, features = thal_genes)
print(p2)
dev.off()

rm(data, data.mouse, data.lamprey, data.mouse.raw, data.lamprey.raw)
