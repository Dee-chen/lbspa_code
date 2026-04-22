# =============================================================================
# 06_multi_species_neuron.R — Multi-species neuronal cell type comparison
# =============================================================================
# Compares neuronal cell types across 8 vertebrate species:
#   Lamprey, Zebrafish, Lizard, Turtle, Pigeon, Zebrafinch, Mouse, Macaque
#
# Workflow:
# 1. Find cell-type markers per species
# 2. Compute cross-species marker intersections and unions
# 3. Generate scaled heatmap with species and cell-type annotations
# 4. Classify genes as species-specific or shared based on expression patterns
# 5. Examine conserved microglia and oligodendrocyte markers
#
# Input:  Multi-species downsampled integrated Seurat object
# Output: Marker CSV files, annotated heatmaps, gene cluster assignments
# =============================================================================

source("00_utils.R")

# --- Load multi-species integrated dataset ------------------------------------
# 8 species, 40000 cells downsampled, 13172 shared orthologs
multispecies <- readRDS(file.path(BASE_DIR,
  "01multispecies_downsample40000_geneIntersect13172_dims20_2000.celltype1119.rds"))

# Standardize species labels
multispecies$species_0 <- as.factor(multispecies$species)
levels(multispecies$species_0) <- c("Lamprey", "Lamprey", "Lizard", "Macaque",
                                     "Mouse", "Pigeon", "Turtle", "Zebrafinch",
                                     "Zebrafish")


# =============================================================================
# SECTION 1: Species-specific marker identification
# =============================================================================
# Find differentially expressed markers for each cell type within each species.

species_list <- c("Lamprey", "Zebrafish", "Lizard", "Turtle",
                  "Pigeon", "Zebrafinch", "Mouse", "Macaque")

markers_by_species <- list()
for (sp in species_list) {
  message("Finding markers for: ", sp)
  sp_data <- subset(multispecies, species_0 == sp)
  markers_by_species[[sp]] <- FindAllMarkers(
    sp_data, group.by = "celltype1119",
    logfc.threshold = 1, min.pct = 0.1)
}

# Save marker lists
for (sp in species_list) {
  write.csv(markers_by_species[[sp]],
            file.path(RESULT_DIR, paste0("00.neurons/marker.", sp, ".csv")))
}


# =============================================================================
# SECTION 2: Neuronal marker union and intersection across species
# =============================================================================
# Focus on neuronal cell types for cross-species comparison

neuron_types <- c("EX", "Granule_like", "IN", "Lizard_PRN",
                  "PRN_1_SNCG", "PRN_2", "Purkinje")

# Filter to upregulated neuronal markers per species
neuron_markers <- lapply(markers_by_species, function(m) {
  subset(m, avg_log2FC > 0 & cluster %in% neuron_types)$gene
})

# Union of all species' neuronal markers
marker_union <- Reduce(union, neuron_markers)
message("Total neuronal marker union: ", length(marker_union), " genes")

# Intersection across all non-Lizard species (7-way)
# (Lizard excluded due to unique PRN population)
non_lizard <- setdiff(species_list, "Lizard")
marker_intersect <- Reduce(intersect, neuron_markers[non_lizard])
message("7-species neuronal marker intersection: ",
        length(marker_intersect), " genes")

# Full 8-species intersection (including all markers, not just neurons)
all_marker_intersect <- Reduce(intersect,
  lapply(markers_by_species, function(m) m$gene))
message("8-species all-marker intersection: ",
        length(all_marker_intersect), " genes")


# =============================================================================
# SECTION 3: Conserved non-neuronal markers (microglia, oligodendrocytes)
# =============================================================================

# --- Microglia conserved markers ---
# Published canonical microglia markers tested across species
published_mic_markers <- c("SPI1", "TARDBP", "IRF8", "CD81", "SALL1", "BIN1",
                           "GRN", "MEF2A", "CST3", "APBB1IP", "MERTK", "HEXB",
                           "P2RY12", "TREM2")

pdf(file.path(RESULT_DIR, "00.neurons/microglia_markers.pdf"),
    width = 12, height = 10)
Stacked_VlnPlot(
  seurat_object = subset(multispecies, celltype1119 == "Microglial"),
  group.by = "species_0",
  features = published_mic_markers, x_lab_rotate = TRUE)
dev.off()

# Find microglia conserved markers using FindConservedMarkers
Idents(multispecies) <- "celltype1119"
mic_conserved <- FindConservedMarkers(
  object = multispecies, ident.1 = "Microglial",
  grouping.var = "species",
  logfc.threshold = 1, min.pct = 0.1)

# --- Oligodendrocyte lineage conserved markers ---
published_oligo_markers <- c("PDGFRA", "CSPG4", "OLIG1", "OLIG2", "SOX10",
                              "MYRF", "ENPP6", "BCAS1", "TCF7L2",
                              "MBP", "PLP1", "MAG", "CNP", "CLDN11", "MAL")

pdf(file.path(RESULT_DIR, "00.neurons/oligo_markers.pdf"),
    width = 14, height = 12)
Stacked_VlnPlot(
  seurat_object = subset(multispecies,
                          celltype1119 %in% c("Oligo", "OPC")),
  group.by = "species_0", split.by = "celltype1119",
  features = published_oligo_markers, x_lab_rotate = TRUE)
dev.off()


# =============================================================================
# SECTION 4: Scaled heatmap with species/celltype annotations
# =============================================================================

# Prepare neuronal subset (exclude Lizard and Mouse-PRN_1_SNCG outlier)
multispecies.neuron <- subset(multispecies,
  celltype1119 %in% neuron_types & species_0 != "Lizard")
multispecies.neuron <- subset(multispecies.neuron,
  species_celltype != "Mouse-PRN_1_SNCG")

# Create combined species-celltype label for grouping
multispecies$species_celltype <- paste(multispecies$species_0,
                                       multispecies$celltype1119, sep = "-")

# Filter to groups with sufficient cells (>= 30)
cell_counts <- table(multispecies$species_celltype)
keep_types <- names(cell_counts[cell_counts >= 30])
multispecies <- subset(multispecies, species_celltype %in% keep_types)

# Compute average expression matrix
avg_exp_data <- AverageExpression(multispecies.neuron,
  features = marker_union, group.by = "species_celltype",
  assay = "RNA")$RNA

# Min-max scale each gene across groups
plot_matrix <- t(apply(avg_exp_data, 1, min_max_scale))

# Build column annotation (species + cell type)
meta_unique <- multispecies.neuron@meta.data %>%
  select(species_celltype, species_0, celltype1119) %>%
  distinct()
meta_unique$match_key <- gsub("_", "-", meta_unique$species_celltype)
rownames(meta_unique) <- meta_unique$match_key
annotation_col <- meta_unique[colnames(plot_matrix),
                               c("species_0", "celltype1119")]

# Annotation colors
species_levels  <- unique(annotation_col$species_0)
celltype_levels <- unique(annotation_col$celltype1119)
cols_species  <- brewer.pal(min(length(species_levels), 8), "Set2")
names(cols_species)  <- species_levels[1:length(cols_species)]
cols_celltype <- brewer.pal(min(length(celltype_levels), 8), "Dark2")
names(cols_celltype) <- celltype_levels[1:length(cols_celltype)]
ann_colors <- list(species_0 = cols_species, celltype1119 = cols_celltype)

# Highlight selected marker genes on heatmap row labels
highlight_genes <- c("SLC17A6", "CACNA2D1", "FOXP1", "CACNA1A", "SST",
                     "GAD2", "NEFM", "TH", "TSHZ2", "NXPH1", "SOX9",
                     "BCL11B", "GRM1", "PTPRK", "SATB2", "LHX1", "TBR1",
                     "GRIK2", "CBLN1", "MYO16")
all_genes <- rownames(plot_matrix)
labels_to_show <- rep("", length(all_genes))
names(labels_to_show) <- all_genes
valid_idx <- match(highlight_genes, all_genes)
valid_idx <- valid_idx[!is.na(valid_idx)]
labels_to_show[valid_idx] <- all_genes[valid_idx]

# Heatmap color
heatmap_color <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100)

pdf(file.path(RESULT_DIR, "00.neurons/multi_species_heatmap.pdf"),
    width = 16, height = 20)
pheatmap(
  avg_exp_data,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  clustering_distance_rows = "correlation",
  clustering_method = "ward.D2",
  cutree_rows = 20,
  scale = "row",
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  color = heatmap_color,
  show_rownames = TRUE,
  labels_row = labels_to_show,
  show_colnames = TRUE,
  border_color = NA,
  fontsize_col = 8,
  treeheight_row = 20,
  treeheight_col = 20,
  legend_breaks = c(-5, 5),
  legend_labels = c("Min", "Max"))
dev.off()


# =============================================================================
# SECTION 5: Gene specificity classification
# =============================================================================
# Classify each marker gene as species-specific or shared based on
# its expression pattern across species.

all_species <- unique(annotation_col$species_0)

#' Classify a gene as species-specific or shared
#'
#' @param gene_expr Named numeric vector (one row of scaled expression matrix)
#' @param col_anno Column annotation data frame with species_0 column
#' @param species_list Vector of species names
#' @return Character: "<Species>_Specific" or "Shared/Ambiguous"
get_gene_category <- function(gene_expr, col_anno, species_list) {
  # Compute max expression per species across all cell types
  species_max <- sapply(species_list, function(sp) {
    cols <- rownames(col_anno)[col_anno$species_0 == sp]
    max(gene_expr[cols])
  })

  high_thresh <- 0.2   # above this: "highly expressed"
  low_thresh  <- 0.01  # below this: "not expressed"

  highly_expressed <- names(species_max)[species_max > high_thresh]
  low_expressed    <- names(species_max)[species_max < low_thresh]

  # Species-specific: expressed in exactly one species, absent in all others
  if (length(highly_expressed) == 1) {
    other <- setdiff(species_list, highly_expressed)
    if (all(other %in% low_expressed)) {
      return(paste0(highly_expressed, "_Specific"))
    }
  }
  return("Shared/Ambiguous")
}

# Classify all genes
gene_categories <- apply(plot_matrix, 1, function(x) {
  get_gene_category(x, annotation_col, all_species)
})

gene_cat_df <- data.frame(
  gene = names(gene_categories),
  category = gene_categories,
  stringsAsFactors = FALSE)

# Separate species-specific and shared genes
specific_genes <- gene_cat_df[gene_cat_df$category != "Shared/Ambiguous", ]
shared_genes   <- gene_cat_df[gene_cat_df$category == "Shared/Ambiguous", ]

# For shared genes, select top 1000 by variance for visualization
shared_variances <- rowVars(plot_matrix[shared_genes$gene, ])
top_shared <- names(sort(shared_variances, decreasing = TRUE))[1:1000]

# Combine for final annotated heatmap
final_gene_df <- rbind(
  specific_genes,
  data.frame(gene = top_shared, category = "Shared_HighVar"))
final_gene_df$category <- factor(final_gene_df$category,
  levels = c(paste0(all_species, "_Specific"), "Shared_HighVar"))
final_gene_df <- final_gene_df[order(final_gene_df$category), ]

# Row annotation for gene specificity blocks
annotation_row <- data.frame(Gene_Group = final_gene_df$category)
rownames(annotation_row) <- final_gene_df$gene

n_groups <- length(unique(annotation_row$Gene_Group))
block_colors <- c(brewer.pal(min(n_groups - 1, 9), "Set1"), "#999999")
names(block_colors) <- levels(final_gene_df$category)
ann_colors$Gene_Group <- block_colors


# =============================================================================
# SECTION 6: Gene cluster extraction from heatmap
# =============================================================================
# Extract hierarchical clustering assignments for downstream analysis

ph <- pheatmap(plot_matrix, cluster_rows = TRUE, cluster_cols = TRUE,
               clustering_distance_rows = "correlation",
               clustering_method = "ward.D2",
               cutree_rows = 20, scale = "row", silent = TRUE)

gene_clusters <- cutree(ph$tree_row, k = 20)
gene_cluster_df <- data.frame(
  Gene = names(gene_clusters),
  Cluster = paste0("Cluster_", gene_clusters))
write.csv(gene_cluster_df,
          file.path(RESULT_DIR, "00.neurons/heatmap_gene_clusters.csv"),
          row.names = FALSE)

# Report visual ordering of clusters (top to bottom in heatmap)
visual_order <- unique(gene_clusters[ph$tree_row$order])
message("Heatmap cluster order (top to bottom): ",
        paste(visual_order, collapse = " -> "))
