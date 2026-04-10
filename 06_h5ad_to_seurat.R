# =============================================================================
# 06_h5ad_to_seurat.R
# =============================================================================
# Step 6 of the pipeline: convert Python scib output h5ads back to Seurat RDS.
#
# Cho et al. — A harmonized single-cell RNA-seq atlas of the human prostate
# Scientific Data (in submission)
#
# Workflow position:
#   04_integration_all_cells.py  → [THIS FILE] → 07_visualization.R
#   05_integration_luminal_cells.py → [THIS FILE] → 07_visualization.R
#
# Why this step is needed:
#   The Python scib pipeline (steps 04 and 05) saves "stripped" h5ad files:
#   count matrix + lognorm layer + all obsm embeddings (scANVI UMAP, etc.),
#   but without the full obs metadata (which would bloat the file and is
#   already saved separately as a CSV).
#   This script re-imports those h5ads via SeuratDisk, re-attaches the full
#   obs metadata from the companion CSV, and saves a complete Seurat RDS
#   that can be used for all downstream R-based plotting (step 07).
#
# How SeuratDisk imports Python h5ads:
#   - Convert(h5ad → h5seurat): translates obsm keys to Seurat reductions,
#     writes lognorm to the RNA "data" slot, raw counts to the "counts" slot.
#   - LoadH5Seurat(): reads the h5seurat into memory as a Seurat object.
#   - obs metadata is NOT included in the stripped h5ad; re-attach from CSV.
#   - UpdateSeuratObject(): upgrades schema to current Seurat version.
#
# Input files (from steps 04 and 05):
#   All-cells (from 04_integration_all_cells.py):
#     merged_a17_hit1_scib_50pc2000g_stripped.h5ad  — embeddings + expression
#     merged_a17_hit1_metadata.csv                  — full obs metadata
#   Luminal (from 05_integration_luminal_cells.py):
#     seu_lum_v1_noNoGroup_scib_30pc3000g_stripped.h5ad  — embeddings + expression
#     seu_lum_v1_noNoGroup_metadata.csv                  — full obs metadata
#
# Output files:
#   merged_a17_hit1_scib_50pc2000g.rds         — all-cells Seurat RDS
#   seu_lum_v1_noNoGroup_scib_30pc3000g.rds    — luminal Seurat RDS
#   (both saved into the scib working directory; used by 07_visualization.R)
#
# Dependencies: Seurat (v5), SeuratDisk
# =============================================================================

library(Seurat)
library(SeuratDisk)

# =============================================================================
# Section 1 — All-cells: merged_a17_hit1_scib_50pc2000g
# =============================================================================
# Source: 04_integration_all_cells.py
#   Integration parameters: 2,000 HVGs, 50 PCs, batch_key = Paper_ID_short,
#   label_key = final_hit1_group_detailed
#   Methods benchmarked: Unintegrated, Harmony, Scanorama, BBKNN, scVI, scANVI
#   Selected method: scANVI (Total scib score = 0.74)
# =============================================================================

setwd("/mctp/share/users/hanbyul/projects/human_prostate_sc/cellranger_7.1/17datasets/python_final")

message("\n=== All-cells: converting h5ad → Seurat RDS ===")

# Step 1: Convert h5ad → h5seurat format (SeuratDisk intermediate)
Convert(
  "merged_a17_hit1_scib_50pc2000g_stripped.h5ad",
  dest     = "h5seurat",
  assay    = "RNA",
  overwrite = TRUE
)

# Step 2: Load h5seurat into Seurat object
#   Only loading the RNA assay (expression data + embeddings).
#   obs metadata columns are empty at this stage; re-attached below.
seu_allcells <- LoadH5Seurat(
  "merged_a17_hit1_scib_50pc2000g_stripped.h5seurat",
  assays = "RNA"
)

message("  Loaded: ", ncol(seu_allcells), " cells, ", nrow(seu_allcells), " genes")
message("  Reductions available: ", paste(names(seu_allcells@reductions), collapse = ", "))

# Step 3: Re-attach full obs metadata from companion CSV
#   The CSV was exported by 04_integration_all_cells.py and contains all
#   original annotation columns (cell type labels, batch info, QC metrics, etc.)
cell_df <- read.csv(
  "merged_a17_hit1_metadata.csv",
  row.names   = 1,
  check.names = FALSE
)
# Align row order to Seurat cell order (in case of any reordering)
cell_df <- cell_df[colnames(seu_allcells), , drop = FALSE]
seu_allcells <- AddMetaData(seu_allcells, metadata = cell_df)

message("  Metadata columns re-attached: ", ncol(cell_df))

# Step 4: Upgrade schema and save RDS
seu_allcells <- UpdateSeuratObject(seu_allcells)
saveRDS(seu_allcells, "merged_a17_hit1_scib_50pc2000g.rds")
message("  RDS saved: merged_a17_hit1_scib_50pc2000g.rds")

rm(seu_allcells); gc()

# =============================================================================
# Section 2 — Luminal cells: seu_lum_v1_noNoGroup_scib_30pc3000g
# =============================================================================
# Source: 05_integration_luminal_cells.py
#   Integration parameters: 3,000 HVGs, 30 PCs, batch_key = Paper_ID_short,
#   label_key = etsgroup
#   Methods benchmarked: Unintegrated, CCA, RPCA, Harmony (from R step 03),
#                        Scanorama, BBKNN, scVI, scANVI (from Python step 05)
#   Selected method: scANVI (Total scib score = 0.78)
#   Luminal ETS subgroups: Benign, BPH, ERG, ETV, nonETS, AR, NEPC, DN
# =============================================================================

setwd("/mctp/share/users/hanbyul/projects/human_prostate_sc/cellranger_7.1/17datasets/python_final/032025_scib")

message("\n=== Luminal cells: converting h5ad → Seurat RDS ===")

# Step 1: Convert h5ad → h5seurat
Convert(
  "seu_lum_v1_noNoGroup_scib_30pc3000g_stripped.h5ad",
  dest      = "h5seurat",
  assay     = "RNA",
  overwrite = TRUE
)

# Step 2: Load h5seurat
seu_lum <- LoadH5Seurat(
  "seu_lum_v1_noNoGroup_scib_30pc3000g_stripped.h5seurat",
  assays = "RNA"
)

message("  Loaded: ", ncol(seu_lum), " cells, ", nrow(seu_lum), " genes")
message("  Reductions available: ", paste(names(seu_lum@reductions), collapse = ", "))

# Step 3: Re-attach full obs metadata from CSV
#   Exported by 05_integration_luminal_cells.py; contains etsgroup, Group,
#   GleasonGrade, Paper_ID_short, and all other annotation columns.
cell_df <- read.csv(
  "seu_lum_v1_noNoGroup_metadata.csv",
  row.names   = 1,
  check.names = FALSE
)
cell_df <- cell_df[colnames(seu_lum), , drop = FALSE]
seu_lum <- AddMetaData(seu_lum, metadata = cell_df)

message("  Metadata columns re-attached: ", ncol(cell_df))

# Step 4: Upgrade schema and save RDS
seu_lum <- UpdateSeuratObject(seu_lum)
saveRDS(seu_lum, "seu_lum_v1_noNoGroup_scib_30pc3000g.rds")
message("  RDS saved: seu_lum_v1_noNoGroup_scib_30pc3000g.rds")

rm(seu_lum); gc()

cat("\nDone. → Next: run 07_visualization.R\n")
cat("Outputs:\n")
cat("  python_final/merged_a17_hit1_scib_50pc2000g.rds\n")
cat("  python_final/032025_scib/seu_lum_v1_noNoGroup_scib_30pc3000g.rds\n")
