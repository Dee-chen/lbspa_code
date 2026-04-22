# =============================================================================
# 05_non_neuron_analysis.R — Non-neuronal cell type analysis
# =============================================================================
# Analyzes non-neuronal cell populations (astroglia, ependymal, choroid plexus)
# in the lamprey brain. Examines conserved astrocyte marker expression and
# investigates regional heterogeneity of glial/ependymal cells along the
# anterior-posterior axis.
#
# Input:  Lamprey 3D brain atlas, non-neuron integrated data, snRNA-seq data
# Output: Marker gene lists, spatial expression plots, violin plots
# =============================================================================

source("00_utils.R")

# --- Load datasets -----------------------------------------------------------

# Non-neuron cross-species integrated data
data <- readRDS(file.path(SCVI_PROC, "15.nonNeuron.rds"))
data <- NormalizeData(data)

# Lamprey 3D brain atlas with updated cluster annotations
data.lamprey.3D <- readRDS(LAMPREY_3D_NEW_PATH)
geneNR <- read.csv(GENE_NR_PATH)

# Cluster metadata for class assignment
clusterdata <- read.csv(file.path(BASE_DIR, "lamprey_color_spatial.csv"))
clusterdata.nonneuron <- subset(clusterdata, class %in% c("Chor", "Epend", "Glia"))

# Assign broad class labels based on cluster annotation
data.lamprey.3D$class <- data.lamprey.3D$anno_whole_cluster.2
levels(data.lamprey.3D$class) <- clusterdata$class

# Subset to non-neuronal populations
data.lamprey.3D.nonneuron <- subset(data.lamprey.3D,
                                     class %in% c("Chor", "Epend", "Glia"))


# --- Published astrocyte marker gene panel ------------------------------------
# Canonical vertebrate astrocyte markers mapped to lamprey gene IDs.
# These genes are used to assess conservation of astrocyte identity in lamprey.
astro_markers <- c(
  "MSTRG.13355",         # SLC1A2 - glutamate transporter
  "MSTRG.13979",         # SLC1A3 - glutamate transporter
  "nbisL1-mrna-12861",   # GLUL   - glutamine synthetase
  "MSTRG.13772",         # SLC38A2 - amino acid transporter
  "nbisL1-mrna-7802",    # SLC38A3 - amino acid transporter
  "MSTRG.21819",         # GFAP/VIM - intermediate filament
  "nbisL1-mrna-13502",   # ALDH1L1 - aldehyde dehydrogenase
  "MSTRG.15566",         # APOB   - apolipoprotein
  "MSTRG.1205",          # FOXG1  - forebrain marker
  "nbisL1-mrna-10174",   # SLC6A11 - GABA transporter
  "MSTRG.30376",         # AQP4   - aquaporin
  "MSTRG.6920"           # HOXC6A - posterior HOX gene
)

# Genes distinguishing forebrain vs hindbrain glial identity
fore_hind_glia_genes <- c(
  "MSTRG.1205",          # FOXG1  - forebrain-enriched
  "MSTRG.2429",
  "nbisL1-mrna-17241",
  "nbisL1-mrna-10174",   # SLC6A11
  "MSTRG.30376",         # AQP4
  "nbisL1-mrna-17385",
  "nbisL1-mrna-17420",
  "MSTRG.6920",          # HOXC6A - hindbrain-enriched
  "MSTRG.7044"
)

# Selected genes for main figure
figure_genes <- c(
  "MSTRG.13355",         # SLC1A2
  "MSTRG.13979",         # SLC1A3
  "nbisL1-mrna-12861",   # GLUL
  "MSTRG.13772",         # SLC38A2
  "nbisL1-mrna-7802",    # SLC38A3
  "MSTRG.1643",          # GPSM2
  "MSTRG.10368",         # NOG2
  "MSTRG.1205",          # FOXG1
  "MSTRG.6920"           # HOXC6A
)


# --- Spatial distribution of glia and ependymal cells -------------------------
# Remove meninges to focus on brain-intrinsic populations
data.lamprey.3D <- subset(data.lamprey.3D, subset = mergeRegion != "Meninges")
Idents(data.lamprey.3D) <- "class"

pdf(file.path(RESULT_DIR, "non_neuron.spatial.pdf"), width = 12, height = 10)

# Spatial maps highlighting glia and ependymal cells (sections 15, 18, 23)
for (sl in c("15", "18", "23")) {
  # Cell class highlights
  p <- Cluster_Highlight_Plot(
    seurat_object = subset(data.lamprey.3D, subset = slices == sl),
    cluster_name = c("Glia", "Epend"), reduction = "align_spatial",
    highlight_color = c("#018E6F", "#6E0000"))
  print(p)

  # Astrocyte marker spatial expression
  p2 <- FeaturePlot_scCustom(
    seurat_object = subset(data.lamprey.3D, slices == sl),
    reduction = "align_spatial", max.cutoff = 5,
    features = figure_genes, num_columns = 3)
  print(p2)
}

# Brain region overview (anterior-posterior axis)
for (sl in c("18", "23")) {
  p <- DimPlot(subset(data.lamprey.3D, subset = slices == sl),
               reduction = "align_spatial", group.by = "level1")
  print(p)
}

dev.off()


# --- Violin plots: astrocyte markers across glial subtypes --------------------

pdf(file.path(RESULT_DIR, "non_neuron.astro_markers.pdf"),
    width = 14, height = 20)

# Astrocyte marker expression across glial clusters
Stacked_VlnPlot(
  seurat_object = subset(data.lamprey.3D, class == "Glia"),
  group.by = "anno_whole_cluster.2",
  features = astro_markers, x_lab_rotate = TRUE)

# Ependymal cell marker expression
Stacked_VlnPlot(
  seurat_object = subset(data.lamprey.3D, class == "Epend"),
  group.by = "anno_whole_cluster.2",
  features = astro_markers, x_lab_rotate = TRUE)

# Forebrain vs hindbrain glial gene panel
Stacked_VlnPlot(
  seurat_object = subset(data.lamprey.3D,
                          class %in% c("Glia", "Epend")),
  group.by = "anno_whole_cluster.2",
  features = fore_hind_glia_genes, x_lab_rotate = TRUE)

dev.off()


# --- Ependymal cell markers --------------------------------------------------
# Identify ependymal-specific markers for the lamprey brain atlas

Idents(data.lamprey.3D) <- "class"
epend.marker <- FindMarkers(data.lamprey.3D, ident.1 = "Epend",
                            logfc.threshold = 0.5, min.pct = 0.1)
epend.marker$geneID <- rownames(epend.marker)

# Annotate with gene names from BLAST NR results
epend.annotated <- merge(x = epend.marker, y = geneNR[, c(4, 6)],
                          by.x = "geneID", by.y = "geneID", all.x = TRUE)
epend.annotated <- subset(epend.annotated, avg_log2FC > 0)
epend.annotated <- epend.annotated[order(-epend.annotated$avg_log2FC), ]
write.csv(epend.annotated, file.path(RESULT_DIR, "10.epend/epend.marker.csv"))
