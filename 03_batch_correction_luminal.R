# =============================================================================
# 03_batch_correction_luminal.R
# =============================================================================
# Step 3 of the pipeline: Seurat v5 batch correction for luminal cells.
# Runs CCA, RPCA, and Harmony integration; exports embeddings + UMAPs as h5ad
# for import into the Python scib benchmarking pipeline.
#
# Cho et al. — A harmonized single-cell RNA-seq atlas of the human prostate
# Scientific Data (in submission)
#
# Workflow position:
#   02_data_merge_and_filter.py → [THIS FILE] → 05_integration_luminal_cells.py
#
# Why this step exists in R (not Python):
#   CCA and RPCA are Seurat-native integration methods (IntegrateLayers).
#   Their corrected embeddings are exported here as h5ad so that
#   05_integration_luminal_cells.py can include them in the scib benchmark
#   alongside the Python-native methods (Scanorama, BBKNN, scVI, scANVI).
#
# Input:
#   merged_a17_STRIPPED_lum_v1_noNoGroup.h5ad — luminal stripped h5ad (from step 02)
#   obs_metadata_lum_v1_noNoGroup.csv          — full obs metadata (from step 02)
#
# Output:
#   seu_lum_v1_noNoGroup.rds   — full Seurat object with all embeddings
#   seu_lum_v1_noNoGroup.h5ad  — h5ad with CCA/RPCA/Harmony embeddings + UMAPs
#                                → imported by 05_integration_luminal_cells.py
#
# Integration methods run here (Seurat v5 IntegrateLayers):
#   1. Unintegrated (PCA baseline, dims 1:30)
#   2. CCA  — Canonical Correlation Analysis
#   3. RPCA — Reciprocal PCA
#   4. Harmony
#   All run on 3,000 HVGs, 50 PCA components (dims 1:30 for kNN and UMAP).
#
# How embeddings appear in the exported h5ad (before Python renaming):
#   X_pca              → X_pca_unintegrated  (unintegrated PCA baseline)
#   X_umap.unintegrated → X_umap_unintegrated
#   X_integrated.cca   → X_cca               (CCA corrected embedding)
#   X_umap.cca         → X_umap_cca
#   X_integrated.rpca  → X_rpca              (RPCA corrected embedding)
#   X_umap.rpca        → X_umap_rpca
#   X_harmony          → X_harmony
#   X_umap.harmony     → X_umap_harmony
# (Renaming is done in 05_integration_luminal_cells.py Step 2)
#
# Dependencies: Seurat (v5), SeuratDisk, SeuratWrappers, reticulate, ggplot2
# =============================================================================

setwd("/mctp/share/users/hanbyul/projects/human_prostate_sc/cellranger_7.1/17datasets/python_final/032025_scib")

library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(SeuratWrappers)
library(reticulate)
library(ggplot2)
library(patchwork)

# Allow large global objects in parallelization
options(future.globals.maxSize = 4e9)

# Use the project Python virtual environment (needed for reticulate / h5ad export)
use_virtualenv(
  "/mctp/share/users/hanbyul/projects/human_prostate_sc/environment_storage/python_venv/human_sc",
  required = TRUE
)

# =============================================================================
# Color palettes (for UMAP diagnostic plots)
# =============================================================================
etsgroup_col <- c(
  "N"="#71A45E","Benign"="#87CEFA","BPH"="#0047ab","nonETS"="#F781BF",
  "ERG"="#984EA3","ETV"="#FF7F00","AR"="#E41A1C","NEPC"="#A65628","DN"="#4DAF4A",
  "No_Group"="#999999"
)
GleasonGrade_col <- setNames(
  c("#71A45E","#87CEFA","#0047ab","#FFEDA0","#FFEDA0","#FEB24C","#FEB24C",
    "#FC4E2A","#E31A1C","#800026","#9E9AC8","#6A51A3","#4A1486"),
  c("N","B","BPH","PCa","GG1","GG2","GG2_GG3","GG3","GG4","GG5","mCRPC","mCRPC-NETD","mCRPC-teSCC")
)
paper_col <- c(
  "P01"="#009CD2","P02"="#e6194B","P03"="#3cb44b","P05"="#ffe119","P06"="#4363d8",
  "P07"="#f58231","P08"="#911eb4","P09"="#46f0f0","P10"="#2b1700","P11"="#f032e6",
  "P12"="#bcf60c","P14"="#fabebe","P15"="#008080","P16"="#e6beff","P17"="#9A6324",
  "P19"="#FFFF33","P21"="#999999"
)

# Batch key (must match Python scripts)
BATCH_KEY  <- "Paper_ID_short"
COLOR_KEYS <- c("etsgroup", "GleasonGrade", "Paper_ID_short")

palette_for_key <- list(
  etsgroup       = etsgroup_col,
  GleasonGrade   = GleasonGrade_col,
  Paper_ID_short = paper_col
)

# =============================================================================
# Helper: load stripped h5ad → Seurat v5
# =============================================================================
load_lum_seurat <- function(h5ad_file, obs_csv, label) {
  message("\n=== Loading ", label, " ===")

  h5seurat_file <- sub("\\.h5ad$", ".h5seurat", h5ad_file)

  # Convert h5ad → h5seurat (SeuratDisk)
  Convert(h5ad_file, dest = "h5seurat", assay = "RNA", overwrite = TRUE)
  seu <- LoadH5Seurat(h5seurat_file)

  # Re-attach full metadata from CSV
  # (obs metadata was stripped from the h5ad in 02_data_merge_and_filter.py
  #  to keep file size manageable for SeuratDisk import)
  cell_df <- read.csv(obs_csv, row.names = 1, check.names = FALSE)
  cell_df  <- cell_df[colnames(seu), , drop = FALSE]
  seu      <- AddMetaData(seu, metadata = cell_df)

  # Upgrade to Seurat v5; convert RNA assay to Assay5 before split()
  seu <- UpdateSeuratObject(seu)
  seu[["RNA"]] <- as(seu[["RNA"]], "Assay5")

  message("  Cells: ", ncol(seu), "  Genes: ", nrow(seu))
  seu
}

# =============================================================================
# Helper: run Seurat v5 integration (CCA, RPCA, Harmony) and export h5ad
# =============================================================================
run_integration <- function(seu, label, out_pdf, out_rds, out_h5ad) {
  message("\n=== Running integration: ", label, " ===")

  # Split RNA assay by batch for Seurat v5 multi-layer integration
  seu[["RNA"]] <- split(seu[["RNA"]], f = seu[[BATCH_KEY, drop = TRUE]])

  # -- Shared preprocessing --------------------------------------------------
  seu <- NormalizeData(seu)
  seu <- FindVariableFeatures(seu, nfeatures = 3000)
  seu <- ScaleData(seu)
  seu <- RunPCA(seu, npcs = 50)

  # -- Unintegrated baseline (PCA, dims 1:30) --------------------------------
  seu <- FindNeighbors(seu, dims = 1:30, reduction = "pca")
  seu <- FindClusters(seu, resolution = 0.5, cluster.name = "unintegrated_clusters")
  seu <- RunUMAP(seu, dims = 1:30, reduction = "pca",
                 reduction.name = "umap.unintegrated")
  message("  Unintegrated: done")

  # -- CCA -------------------------------------------------------------------
  # Canonical Correlation Analysis: finds shared correlation structure across batches
  seu <- IntegrateLayers(seu, method = CCAIntegration,
                         orig.reduction = "pca",
                         new.reduction  = "integrated.cca",
                         verbose = FALSE)
  seu <- FindNeighbors(seu, reduction = "integrated.cca", dims = 1:30)
  seu <- FindClusters(seu, resolution = 0.5, cluster.name = "cca_clusters")
  seu <- RunUMAP(seu, reduction = "integrated.cca", dims = 1:30,
                 reduction.name = "umap.cca")
  message("  CCA: done")

  # -- RPCA ------------------------------------------------------------------
  # Reciprocal PCA: projects each dataset into the PCA space of the other;
  # more conservative than CCA (better preserves within-batch structure)
  seu <- IntegrateLayers(seu, method = RPCAIntegration,
                         orig.reduction = "pca",
                         new.reduction  = "integrated.rpca",
                         verbose = FALSE)
  seu <- FindNeighbors(seu, reduction = "integrated.rpca", dims = 1:30)
  seu <- FindClusters(seu, resolution = 0.5, cluster.name = "rpca_clusters")
  seu <- RunUMAP(seu, reduction = "integrated.rpca", dims = 1:30,
                 reduction.name = "umap.rpca")
  message("  RPCA: done")

  # -- Harmony ---------------------------------------------------------------
  seu <- IntegrateLayers(seu, method = HarmonyIntegration,
                         orig.reduction = "pca",
                         new.reduction  = "harmony",
                         verbose = FALSE)
  seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:30)
  seu <- FindClusters(seu, resolution = 0.5, cluster.name = "harmony_clusters")
  seu <- RunUMAP(seu, reduction = "harmony", dims = 1:30,
                 reduction.name = "umap.harmony")
  message("  Harmony: done")

  # -- Diagnostic DimPlots → PDF --------------------------------------------
  methods_info <- list(
    list(reduction = "umap.unintegrated", label = "Unintegrated"),
    list(reduction = "umap.cca",          label = "CCA"),
    list(reduction = "umap.rpca",         label = "RPCA"),
    list(reduction = "umap.harmony",      label = "Harmony")
  )

  pdf(out_pdf, width = 5 * length(COLOR_KEYS), height = 5)
  for (m in methods_info) {
    plots <- lapply(COLOR_KEYS, function(key) {
      vals <- as.character(unique(seu[[key, drop = TRUE]]))
      cols <- palette_for_key[[key]]
      cols <- cols[names(cols) %in% vals]
      p <- DimPlot(seu, reduction = m$reduction, group.by = key,
                   label = FALSE, raster = FALSE) +
        ggtitle(paste0(m$label, " — ", key)) +
        theme(legend.text = element_text(size = 7),
              legend.key.size = unit(0.4, "cm"))
      if (length(cols) > 0)
        p <- p + scale_color_manual(values = cols, na.value = "lightgrey")
      p
    })
    print(wrap_plots(plots, nrow = 1))
  }
  dev.off()
  message("  PDF saved: ", out_pdf)

  # -- Join layers → save RDS -----------------------------------------------
  seu <- JoinLayers(seu)
  saveRDS(seu, out_rds)
  message("  RDS saved: ", out_rds)

  # -- Export h5ad for Python scib benchmarking ------------------------------
  # DietSeurat drops scale.data (too large for SaveH5Seurat on 130k+ cells).
  # Only the embeddings needed by 05_integration_luminal_cells.py are kept.
  #
  # Embeddings exported (Python will rename them in Step 2 of script 05):
  #   pca              → X_pca_unintegrated
  #   umap.unintegrated → X_umap_unintegrated
  #   integrated.cca   → X_cca
  #   umap.cca         → X_umap_cca
  #   integrated.rpca  → X_rpca
  #   umap.rpca        → X_umap_rpca
  #   harmony          → X_harmony
  #   umap.harmony     → X_umap_harmony
  seu_h5 <- DietSeurat(
    seu,
    layers    = c("counts", "data"),
    dimreducs = c("pca",
                  "umap.unintegrated",
                  "integrated.cca",  "umap.cca",
                  "integrated.rpca", "umap.rpca",
                  "harmony",         "umap.harmony")
  )
  # Downgrade Assay5 → Assay (SeuratDisk only supports Assay)
  seu_h5[["RNA"]] <- as(seu_h5[["RNA"]], "Assay")
  h5seurat_out <- sub("\\.h5ad$", ".h5seurat", out_h5ad)
  SaveH5Seurat(seu_h5, filename = h5seurat_out, overwrite = TRUE)
  Convert(h5seurat_out, dest = "h5ad", overwrite = TRUE)
  message("  h5ad saved: ", out_h5ad)
  message("  → Next: run 05_integration_luminal_cells.py")

  invisible(seu)
}

# =============================================================================
# Run integration for luminal cells v1 (No_Group removed)
# =============================================================================
seu_v1 <- load_lum_seurat(
  h5ad_file = "merged_a17_STRIPPED_lum_v1_noNoGroup.h5ad",
  obs_csv   = "obs_metadata_lum_v1_noNoGroup.csv",
  label     = "v1_noNoGroup"
)

run_integration(
  seu      = seu_v1,
  label    = "v1_noNoGroup",
  out_pdf  = "umap_batch_lum_v1_noNoGroup.pdf",
  out_rds  = "seu_lum_v1_noNoGroup.rds",
  out_h5ad = "seu_lum_v1_noNoGroup.h5ad"
)

rm(seu_v1); gc()
cat("\nDone. → Next: run 05_integration_luminal_cells.py\n")
