# =============================================================================
# 02_conserved_marker_analysis.R — Cross-species conserved marker identification
# =============================================================================
# For each brain region, identifies conserved marker genes between lamprey and
# mouse. Splits integrated data by species, maps fine cell types to broader
# classes, finds markers per class, and computes cross-species marker
# intersections.
#
# Input:  Filtered Seurat objects from 01_region_integration_analysis.R
# Output: Conserved marker CSV files; all-marker lists per region
# =============================================================================

source("00_utils.R")

# =============================================================================
# SECTION 1: Conserved markers at cluster level (integrated clusters)
# =============================================================================
# Identifies markers that are conserved between species within each integrated
# cluster. Uses FindConservedMarkers with species as the grouping variable.

# --- Region-specific cluster configurations ---
conserved_marker_config <- list(
  list(
    name     = "olfactory_bulb",
    input    = file.path(RESULT_DIR, "01.olfactory_bulb.rds"),
    output   = file.path(RESULT_DIR, "01.olfactory_bulb.conservedmarker.csv"),
    clusters = c(0:8, 22)
  ),
  list(
    name     = "cortex",
    input    = file.path(RESULT_DIR, "02.cortex.rds"),
    output   = file.path(RESULT_DIR, "02.cortex.conservedmarker.csv"),
    clusters = c(0:8, 22)
  ),
  list(
    name     = "hippocampal",
    input    = file.path(RESULT_DIR, "03.hippocampal.rds"),
    output   = file.path(RESULT_DIR, "03.hippocampal.conservedmarker.csv"),
    clusters = c(0:8, 52)
  ),
  list(
    name     = "striatum",
    input    = file.path(RESULT_DIR, "04.striatum.rds"),
    output   = file.path(RESULT_DIR, "04.striatum.conservedmarker.csv"),
    clusters = 0:2
  ),
  list(
    name     = "hypothalamus",
    input    = file.path(RESULT_DIR, "05.hypothalamus.rds"),
    output   = file.path(RESULT_DIR, "05.hypothalamus.conservedmarker.csv"),
    clusters = c(0:6, 8)
  ),
  list(
    name     = "thalamus",
    input    = file.path(RESULT_DIR, "06.thalamus.rds"),
    output   = file.path(RESULT_DIR, "06.thalamus.conservedmarker.csv"),
    clusters = 0:5
  ),
  list(
    name     = "midbrain",
    input    = file.path(RESULT_DIR, "07.midbrain.rds"),
    output   = file.path(RESULT_DIR, "07.midbrain.conservedmarker.csv"),
    clusters = 0:9
  ),
  list(
    name     = "pons",
    input    = file.path(RESULT_DIR, "08.pons.rds"),
    output   = file.path(RESULT_DIR, "08.pons.conservedmarker.csv"),
    clusters = c(0:11, 56)
  )
)

for (cfg in conserved_marker_config) {
  message("=== Conserved markers: ", cfg$name, " ===")

  data <- readRDS(cfg$input)
  Idents(data) <- "RNA_snn_res.0.5"

  # Find conserved markers across species for each cluster
  marker <- run_conserved_marker_loop(data, cfg$clusters)
  write.csv(marker, cfg$output)
  message("Saved: ", cfg$output)

  # Also find all markers (not species-stratified) for reference
  all_marker <- FindAllMarkers(data, logfc.threshold = 0.1, min.pct = 0.25)
  write.csv(all_marker, gsub("conservedmarker", "allmarker", cfg$output))
}


# =============================================================================
# SECTION 2: Class-level cross-species marker analysis
# =============================================================================
# Splits data into mouse/lamprey, maps fine cell types to broad classes,
# finds class-level markers per species, then intersects them to identify
# evolutionarily conserved cell-type markers.

# --- Cell type to class mapping tables per region ---
# Each region has a unique mapping reflecting the cell composition in that area.
# species==1: mouse; species==0: lamprey

class_mapping_config <- list(

  # --- Olfactory bulb ---
  list(
    name = "olfactory_bulb",
    input = file.path(SCVI_PROC, "01.olfactory_bulb.rds"),
    # Mouse celltype indices and class mapping
    mouse_ct_idx   = c(1, 6, 7, 34:37, 67:69),
    mouse_classes  = c("Astro", "Endo", "Epend", "Exc", "Inh", "OPC", "Oligo",
                       "Exc", "Inh", "Endo"),
    # Lamprey celltype indices and class mapping
    lamprey_ct_idx   = 16:33,
    lamprey_classes  = c("Glia-like", "Glia-like", "Glia-like", "Inh-like",
                         "Glia-like", "Glia-like", "Glia-like", "Inh-like",
                         "Exc-like", "Inh-like", "Inh-like", "Exc-like",
                         "Inh-like", "Exc-like", "Glia-like", "Glia-like",
                         "Inh-like", "Glia-like")
  ),

  # --- Cortex ---
  list(
    name = "cortex",
    input = file.path(SCVI_PROC, "02.cortex.rds"),
    mouse_ct_idx   = c(1:11, 31:33),
    mouse_classes  = c("Astro", "Endo", "Gran-DG", "Endo", "Epend", "Epend",
                       "Mic", "Exc", "Inh", "OPC", "Oligo", "Exc", "Inh", "Endo"),
    lamprey_ct_idx   = 12:30,
    lamprey_classes  = c("Astro-like", "Exc-like2", "Exc-like2", "Astro-like",
                         "Exc-like1", "Exc-like3", "Exc-like3", "Exc-like3",
                         "Epend-like", "Astro-like", "Exc-like2", "Exc-like3",
                         "Exc-like1", "Exc-like2", "Inh-like", "Inh-like",
                         "Inh-like", "Astro-like", "Epend-like")
  ),

  # --- Hippocampal ---
  list(
    name = "hippocampal",
    input = file.path(SCVI_PROC, "03.Hippocampal.rds"),
    mouse_ct_idx   = c(1, 3, 5, 6, 8, 11, 12, 27:29),
    mouse_classes  = c("Astro", "Gran", "Endo", "Epend", "Micro", "OPC",
                       "Oligo", "Exc", "Inh", "Endo"),
    lamprey_ct_idx   = 13:26,
    lamprey_classes  = c("Glia-like", "Glia-like", "Glia-like", "Glia-like",
                         "Exc-like", "Exc-like", "Inh-like", "Inh-like",
                         "Inh-like", "Exc-like", "Exc-like", "Exc-like",
                         "Inh-like", "Inh-like")
  ),

  # --- Striatum ---
  list(
    name = "striatum",
    input = file.path(SCVI_PROC, "04.Striatum.rds"),
    mouse_ct_idx   = c(1, 7, 8, 10:14, 23:25),
    mouse_classes  = c("Astro", "Gran", "Endo", "Epend", "Micro", "OPC",
                       "Oligo", "Exc", "Inh", "Endo"),
    lamprey_ct_idx   = 15:22,
    lamprey_classes  = c("Glia-like", "Glia-like", "Glia-like", "Glia-like",
                         "Exc-like", "Exc-like", "Inh-like", "Inh-like",
                         "Inh-like", "Exc-like", "Exc-like", "Exc-like",
                         "Inh-like", "Inh-like")
  ),

  # --- Hypothalamus ---
  list(
    name = "hypothalamus",
    input = file.path(SCVI_PROC, "05.Hypothalamus.rds"),
    mouse_ct_idx   = c(1:6, 27:30, 33),
    mouse_classes  = c("Astro", "Mono", "Chor", "DiMes", "Endo", "Epend",
                       "Epend", "Micro", "OPC", "Oligo", "Endo"),
    lamprey_ct_idx   = 7:26,
    lamprey_classes  = c("Glia-like", "Glia-like", "Glia-like", "DiMes-like",
                         "DiMes-like", "Mono-like", "DiMes-like", "Mono-like",
                         "DiMes-like", "DiMes-like", "HypNeurons", "Glia-like",
                         "HypNeurons", "Mono-like", "DiMes-like", "Epend-like",
                         "Glia-like", "Glia-like", "Epend-like", "Epend-like")
  ),

  # --- Thalamus ---
  list(
    name = "thalamus",
    input = file.path(SCVI_PROC, "06.Thalamus.rds"),
    mouse_ct_idx   = c(1:8, 10, 11, 32),
    mouse_classes  = c("Astro", "Mono", "Chor", "DiMes", "Endo", "Epend",
                       "Epend", "Micro", "OPC", "Oligo", "Endo"),
    lamprey_ct_idx   = 15:31,
    lamprey_classes  = c("Glia-like", "Glia-like", "DiMes-like", "Mono-like",
                         "Mono-like", "DiMes-like", "Endo-like", "Endo-like",
                         "Endo-like", "Mono-like", "DiMes-like", "DiMes-like",
                         "Mono-like", "Mono-like", "DiMes-like", "Mono-like",
                         "Glia-like")
  ),

  # --- Midbrain ---
  list(
    name = "midbrain",
    input = file.path(SCVI_PROC, "07.Midbrain.rds"),
    mouse_ct_idx   = c(1:3, 5:9, 32:34, 36),
    mouse_classes  = c("Astro", "Mono", "Chor", "DiMes", "Endo", "Epend",
                       "Epend", "Mic", "OPC", "Oligo", "Rh", "Endo"),
    lamprey_ct_idx   = 10:30,
    lamprey_classes  = c("Glia-like1", "MidNeuron", "MidNeuron", "MidNeuron",
                         "DiMes-like1", "Mono-like", "Rh-like", "MidNeuron",
                         "Rh-like", "DiMes-like1", "DiMes-like1", "Mono-like",
                         "Rh-like", "Mono-like", "Mono-like", "DiMes-like1",
                         "MidNeuron", "MidNeuron", "Glia-like2", "Glia-like1",
                         "Glia-like2")
  ),

  # --- Pons ---
  list(
    name = "pons",
    input = file.path(SCVI_PROC, "08.Pons.rds"),
    mouse_ct_idx   = c(1, 17:28, 50:53),
    mouse_classes  = c("Astro", "Neuroblast", "Mono", "Chor", "DiMes", "Endo",
                       "Epend", "Epend", "Mic", "OPC", "Oligo", "Rh",
                       "Exc/Inh", "Exc/Inh", "Endo"),
    lamprey_ct_idx   = 2:16,
    lamprey_classes  = c("Astro-like", "Epend-like", "Epend-like", "Epend-like",
                         "Exc/Inh-like", "Neuroblast", "Rh-like", "Neuroblast",
                         "Neuroblast", "Neuroblast", "Mono-like", "Rh-like",
                         "Mono-like", "Mono-like", "Exc/Inh-like", "Rh-like",
                         "Epend-like", "Astro-like", "Exc/Inh-like", "Mono-like")
  )
)


# --- Run class-level marker analysis for each region -------------------------
for (cfg in class_mapping_config) {
  message("=== Class-level markers: ", cfg$name, " ===")

  data <- readRDS(cfg$input)
  data <- NormalizeData(data)

  # Split by species
  ct_names <- names(table(data$celltype))
  data.mouse   <- subset(data, subset = species == 1 &
                           celltype %in% ct_names[cfg$mouse_ct_idx])
  data.lamprey <- subset(data, subset = species == 0 &
                           celltype %in% ct_names[cfg$lamprey_ct_idx])

  # Map fine cell types to broad classes
  data.mouse$class   <- data.mouse$celltype
  levels(data.mouse$class)   <- cfg$mouse_classes
  data.lamprey$class <- data.lamprey$celltype
  levels(data.lamprey$class) <- cfg$lamprey_classes

  # Find markers per class within each species
  marker.mouse <- FindAllMarkers(data.mouse, group.by = "class",
                                 logfc.threshold = 0.25, min.pct = 0.1)
  marker.lamprey <- FindAllMarkers(data.lamprey, group.by = "class",
                                   logfc.threshold = 0.25, min.pct = 0.1)

  # Filter for upregulated markers
  marker.mouse   <- subset(marker.mouse, avg_log2FC > 0)
  marker.lamprey <- subset(marker.lamprey, avg_log2FC > 0)

  # Identify cross-species conserved markers (gene name intersection)
  marker.union <- intersect(marker.mouse$gene, marker.lamprey$gene)
  message("  Conserved marker genes: ", length(marker.union))

  # Generate heatmaps for mouse and lamprey
  output_prefix <- file.path(RESULT_DIR, cfg$name)

  pdf(paste0(output_prefix, ".class_heatmap.pdf"), width = 10, height = 8)
  plot_heatmap_scaled(data.mouse, marker.union, "class",
                      color_palette = "#C1640E",
                      title = paste(cfg$name, "- Mouse"))
  plot_heatmap_scaled(data.lamprey, marker.union, "class",
                      color_palette = "#270A60",
                      title = paste(cfg$name, "- Lamprey"))
  dev.off()

  # Save marker intersection
  write.csv(data.frame(gene = marker.union),
            paste0(output_prefix, ".conserved_class_markers.csv"),
            row.names = FALSE)
}


# =============================================================================
# SECTION 3: Midbrain-specific curated conserved marker list
# =============================================================================
# Hand-curated conserved markers for the midbrain, used for targeted heatmaps.

midbrain_conserved_genes <- c(
  # Glial markers
  "ALDOC", "SLC1A2", "SPARC", "CPE",
  # Pan-neuronal markers
  "RTN1", "STMN2", "VAMP2", "HSPE1", "ELAVL2",
  # Excitatory / glutamatergic markers
  "PSMD12", "CPLX1", "PDHB",
  "CNTNAP2", "SNAP25", "SLC17A6", "GRIA2", "MDH1", "VSNL1", "PPP3R1", "KCNC1",
  # Microglia / immune-related
  "TUBB4B", "CRIP2", "LGMN", "LRP1",
  # Oligodendrocyte lineage
  "CNP", "ERMN", "ATP1B3", "PHGDH", "ANK3",
  # Additional neuronal markers
  "NEFL", "PRDX5", "GOT1", "PKM", "CYCS", "SLC25A3"
)

# This list is also reused for the pons analysis (same conserved gene set)
pons_conserved_genes <- c(midbrain_conserved_genes, "NRXN1", "CHCHD6", "MEA1")
