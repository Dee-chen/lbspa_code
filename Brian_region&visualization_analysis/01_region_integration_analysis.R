# =============================================================================
# 01_region_integration_analysis.R — scVI integration and Seurat preprocessing
# =============================================================================
# Processes scVI-corrected lamprey-mouse integrated data for each brain region.
# For each region: reads integrated object -> standard Seurat pipeline ->
# quality filtering -> saves cleaned object.
#
# Input:  scVI-corrected RDS files (one per brain region)
# Output: Filtered Seurat objects saved to RESULT_DIR
# =============================================================================

source("00_utils.R")

# --- Region configuration ----------------------------------------------------
# Each entry: region name, input path, output filename, clusters to keep
region_config <- list(
  list(
    name     = "olfactory_bulb",
    input    = file.path(SCVI_RAW, "01.olfactory_bulb/scvi_corrected.rds"),
    output   = file.path(RESULT_DIR, "01.olfactory_bulb.rds"),
    keep_clusters = c(0:8, 22)
  ),
  list(
    name     = "cortex",
    input    = file.path(SCVI_RAW, "02.cortex/scvi_corrected_1.rds"),
    output   = file.path(RESULT_DIR, "02.cortex.rds"),
    keep_clusters = c(0:9, 52)
  ),
  list(
    name     = "hippocampal",
    input    = file.path(SCVI_RAW, "03.hippocampal/scvi_corrected_2.rds"),
    output   = file.path(RESULT_DIR, "03.hippocampal.rds"),
    keep_clusters = c(0:9, 52)
  ),
  list(
    name     = "striatum",
    input    = file.path(SCVI_RAW, "04.striatum/scvi_corrected_3.rds"),
    output   = file.path(RESULT_DIR, "04.striatum.rds"),
    keep_clusters = c(0:9, 52)
  ),
  list(
    name     = "hypothalamus",
    input    = file.path(SCVI_RAW, "05.hypothalamus/scvi_corrected_4.rds"),
    output   = file.path(RESULT_DIR, "05.hypothalamus.rds"),
    keep_clusters = c(0:6, 8)
  ),
  list(
    name     = "thalamus",
    input    = file.path(SCVI_RAW, "06.thalamus/scvi_corrected_5.rds"),
    output   = file.path(RESULT_DIR, "06.thalamus.rds"),
    keep_clusters = 0:5
  ),
  list(
    name     = "midbrain",
    input    = file.path(SCVI_RAW, "07.midbrain/scvi_corrected_6.rds"),
    output   = file.path(RESULT_DIR, "07.midbrain.rds"),
    keep_clusters = 0:9
  ),
  list(
    name     = "pons",
    input    = file.path(SCVI_RAW, "08.pons/scvi_corrected_7.rds"),
    output   = file.path(RESULT_DIR, "08.pons.rds"),
    keep_clusters = c(0:11, 56)
  )
)

# --- Run integration pipeline for each region --------------------------------
for (cfg in region_config) {
  message("=== Processing: ", cfg$name, " ===")

  # Read scVI-corrected integrated data
  data <- readRDS(cfg$input)

  # Standard Seurat preprocessing on scVI-aligned embeddings
  data <- preprocess_scvi_data(data)

  # Quality check: cell distribution across species per cluster
  message("Cluster x species distribution:")
  print(table(data$RNA_snn_res.0.5, data$species))

  # UMAP overview (one representative plot per region)
  pdf(file.path(RESULT_DIR, paste0(cfg$name, ".umap_overview.pdf")),
      width = 10, height = 5)
  p <- DimPlot(data, group.by = "RNA_snn_res.0.5", label = TRUE,
               split.by = "species")
  print(p)
  dev.off()

  # Filter to retain high-quality clusters
  data <- subset(data, subset = RNA_snn_res.0.5 %in% cfg$keep_clusters)

  # Save filtered object
  saveRDS(data, cfg$output)
  message("Saved: ", cfg$output, " (", ncol(data), " cells)\n")
}
