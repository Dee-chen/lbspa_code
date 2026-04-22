# =============================================================================
# 04_pineal_cross_species.R — Pineal gland cross-species comparison
# =============================================================================
# Compares pineal cell types across three species: lamprey (spatial
# transcriptomics), zebrafish (scRNA-seq), and rat (scRNA-seq).
# Identifies conserved pineal marker genes and examines melatonin
# biosynthesis pathway genes (AANAT, ASMT).
#
# Input:  Lamprey 3D spatial object, zebrafish & rat pineal scRNA-seq
# Output: Cross-species marker comparisons, expression plots
# =============================================================================

source("00_utils.R")

# --- Load datasets -----------------------------------------------------------

# Lamprey pineal spatial data (stereo-seq)
data.lamprey <- readRDS(file.path(BASE_DIR, "17_pin/lamrepy_brain_stereo_3d.rds"))

# Lamprey 3D brain atlas — subset to pineal region
data.lamprey.3D <- readRDS(LAMPREY_3D_NEW_PATH)
data.lamprey.3D.pin <- subset(data.lamprey.3D,
                               subset = mergeRegion == "Pineal")

# Zebrafish and rat pineal scRNA-seq datasets
fish.pineal <- readRDS(file.path(BASE_DIR, "17_pin/fish_brain.rds"))
rat.pineal  <- readRDS(file.path(BASE_DIR, "17_pin/rat_brain.rds"))

# Gene annotation (lamprey gene ID to ortholog name mapping)
geneNR <- read.csv(GENE_NR_PATH)
lamprey.gene <- cbind(rownames(data.lamprey),
                      data.lamprey@assays$RNA@meta.features)

# --- Normalize all datasets ---------------------------------------------------
fish.pineal        <- NormalizeData(fish.pineal)
rat.pineal         <- NormalizeData(rat.pineal)
data.lamprey.3D.pin <- NormalizeData(data.lamprey.3D.pin)


# --- Find species-specific markers -------------------------------------------

# Lamprey pineal cluster markers
lamprey.marker <- FindAllMarkers(data.lamprey.3D.pin,
                                  group.by = "anno_whole_cluster.2",
                                  logfc.threshold = 0.25, min.pct = 0.1)
lamprey.marker <- subset(lamprey.marker, avg_log2FC > 0)

# Map lamprey gene IDs to ortholog names for cross-species comparison
lamprey.marker.name <- subset(lamprey.gene,
                               ids %in% lamprey.marker$gene)

# Zebrafish and rat pineal markers
fish.marker <- FindAllMarkers(fish.pineal, group.by = "celltype",
                               logfc.threshold = 1, min.pct = 0.1)
fish.marker <- subset(fish.marker, avg_log2FC > 1)

rat.marker <- FindAllMarkers(rat.pineal, group.by = "celltype",
                              logfc.threshold = 1, min.pct = 0.1)
rat.marker <- subset(rat.marker, avg_log2FC > 1)


# --- Three-way cross-species marker intersection -----------------------------
# Step 1: lamprey-rat intersection
union.lamprey.rat <- intersect(
  lamprey.marker.name$`rownames(data.lamprey)`,
  rat.marker$gene)

# Step 2: further intersect with zebrafish
marker.three_species <- intersect(union.lamprey.rat, fish.marker$gene)
message("Three-species conserved pineal markers: ",
        length(marker.three_species))


# --- Photoreceptor-specific comparison: lamprey Photo-3/5 vs rat b-Pinealocyte ---
# Lamprey pineal photoreceptor clusters 3 and 5
lamprey.photo.gene <- subset(
  lamprey.gene,
  ids %in% subset(lamprey.marker,
                  cluster %in% c("Pin-Photo-3", "Pin-Photo-5"))$gene)

# Intersect with rat b-Pinealocyte markers
photo_vs_bpin <- intersect(
  lamprey.photo.gene$`rownames(data.lamprey)`,
  subset(rat.marker, cluster == "b-Pinealocyte")$gene)
photo_vs_bpin.ids <- lamprey.gene[photo_vs_bpin, ]
message("Photo-3/5 vs b-Pinealocyte shared markers: ", length(photo_vs_bpin))


# --- Melatonin biosynthesis pathway genes -------------------------------------
# AANAT and ASMT are key enzymes in melatonin synthesis.
# Check their expression across all three species.

# Lamprey orthologs of AANAT and ASMT
lamprey_melatonin_genes <- c("nbisL1-mrna-340", "MSTRG.3007")  # ASMT, AANAT
mammal_melatonin_genes  <- c("ASMT", "AANAT")

pdf(file.path(RESULT_DIR, "pineal.melatonin_genes.pdf"),
    width = 10, height = 6)

# Rat pineal
p1 <- FeaturePlot_scCustom(seurat_object = rat.pineal,
                            reduction = "umap", max.cutoff = 5,
                            features = mammal_melatonin_genes)
print(p1)

# Lamprey pineal (spatial section 19)
p2 <- FeaturePlot_scCustom(
  seurat_object = subset(data.lamprey.3D.pin, slices == "19"),
  reduction = "align_spatial", max.cutoff = 5,
  features = lamprey_melatonin_genes)
print(p2)

# Zebrafish pineal
p3 <- FeaturePlot_scCustom(seurat_object = fish.pineal,
                            reduction = "umap", max.cutoff = 5,
                            features = mammal_melatonin_genes)
print(p3)

dev.off()


# --- Heatmap: three-species conserved markers ---------------------------------
if (length(marker.three_species) > 0) {
  pdf(file.path(RESULT_DIR, "pineal.conserved_heatmap.pdf"),
      width = 10, height = 8)
  plot_heatmap_scaled(rat.pineal, marker.three_species, "celltype",
                      color_palette = "#C1640E",
                      title = "Pineal conserved markers - Rat")
  dev.off()
}
