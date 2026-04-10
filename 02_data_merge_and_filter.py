"""
02_data_merge_and_filter.py
============================
Step 2 of the pipeline: load per-sample 10x h5 files into Python, merge
into a single AnnData object, filter cells using the R annotation table, and
create the input files for the batch correction steps.

Cho et al. — A harmonized single-cell RNA-seq atlas of the human prostate
Scientific Data (in submission)

Workflow position:
  01_QC_and_annotation.R → [THIS FILE] → 03_batch_correction_luminal.R
                                        → 04_integration_all_cells.py

What this script does:
  Part A — Per-sample loading and QC metric export (pre-filter):
    1. Load each Cell Ranger v7.1 filtered_feature_bc_matrix.h5 per sample
       (P02/Wong 2022 handled separately using HTO demultiplexing barcodes)
    2. Compute per-cell QC metrics (nCount, nFeature, pct_counts_mt)
    3. Save pre-filter QC obs table → QC_metrics_raw_prefilter.csv
       (used by 05_visualization.R for Supplementary Figure 1a)
    4. Save per-study h5ad files → {study}_combined.h5ad
    5. Concatenate all studies → all_studies_combined.h5ad

  Part B — Filter by R cell type annotation:
    6. Load all_studies_combined.h5ad + metacomb_17data_final.csv (from R step 1)
    7. Keep only cells present in the annotation table; attach R annotations
    8. Save → merged_a17.h5ad   (all annotated cells, raw counts, all genes)

  Part C — Create integration inputs:
    9. Remove low-quality / unannotated cells; remove MT genes
       → merged_a17_hit1.h5ad   (input for 04_integration_all_cells.py)
    10. Subset to epithelial cells only; strip to bare h5ad (no scale.data)
        → merged_a17_STRIPPED_lum_v1_noNoGroup.h5ad
        → obs_metadata_lum_v1_noNoGroup.csv
        (inputs for 03_batch_correction_luminal.R)

Key output files:
  python/data/{study}_combined.h5ad         — per-study raw h5ad
  python_final/QC_metrics_raw_prefilter.csv — pre-filter QC metrics (for R plotting)
  python_final/merged_a17.h5ad              — all cells, post annotation merge
  python_final/merged_a17_hit1.h5ad         — filtered cells, input for all-cells integration
  python_final/032025_scib/merged_a17_STRIPPED_lum_v1_noNoGroup.h5ad — luminal input for R
  python_final/032025_scib/obs_metadata_lum_v1_noNoGroup.csv

Dependencies: scanpy, anndata, pandas, numpy, scipy
"""

import os
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc

sc.settings.verbosity = 3

# =============================================================================
# Paths
# =============================================================================
BASE_DIR   = "/mctp/share/users/hanbyul/projects/human_prostate_sc/cellranger_7.1/17datasets"
DATA_DIR   = os.path.join(BASE_DIR, "python/data")          # per-study h5ad cache
PY_FINAL   = os.path.join(BASE_DIR, "python_final")         # merged outputs
SCIB_DIR   = os.path.join(PY_FINAL, "032025_scib")          # luminal outputs for R

os.makedirs(DATA_DIR, exist_ok=True)
os.makedirs(PY_FINAL, exist_ok=True)
os.makedirs(SCIB_DIR, exist_ok=True)

# Study order (consistent with R paperorder)
study_ids = [
    "P16_2018_CellRep_Henry",  "P17_2022_JPathol_Joseph",  "P19_2024_BioRxiv_TSC",
    "P14_2020_Elife_Crowley",  "P15_2021_JPathol_Joseph",  "P03_2020_MolCancer_Ma",
    "P01_2021_CellRep_Tuong",  "P05_2023_NatCommun_Hirz",  "P02_2022_NatCommun_Wong",
    "P10_2024_NatCommun_Graham","P06_2022_MolCancer_Heidegger","P11_2022_JExpMed_Masetti",
    "P21_2024_NPJPrecisOncol_Cheng","P08_2020_CommunBiol_Dong","P09_2022_EurUrol_Cheng",
    "P12_2021_CancerCell_Kfoury","P07_2022_Science_Chan",
]

# Cells to exclude in Part C (low-quality / unannotated labels)
EXCLUDE_LABELS = {"lowqual", "unknown", "remove"}

# =============================================================================
# Part A — Per-sample loading, QC metric export, per-study h5ad files
# =============================================================================
print("\n" + "="*60)
print("Part A: Load per-sample h5 files, compute QC, save per-study h5ad")
print("="*60)

meta_clean = pd.read_csv(os.path.join(BASE_DIR, "meta.clean.csv"), index_col=0)

# Wong 2022: load HTO barcode assignments exported from R
# (wongbarcode.csv exported from wong.tumor.benign.barcodes.rds in R)
wongbarcode = pd.read_csv(os.path.join(BASE_DIR, "wongbarcode.csv"))
wb_by_patient = {pid: grp.set_index("barcode") for pid, grp in wongbarcode.groupby("patient_id")}

all_qc_dfs = []  # collect pre-filter QC obs across all studies

for study in study_ids:
    out_path = os.path.join(DATA_DIR, f"{study}_combined.h5ad")
    if os.path.exists(out_path):
        print(f"[skip] {study} already processed")
        # Still collect QC obs for the CSV export
        a_tmp = sc.read_h5ad(out_path, backed="r")
        obs_tmp = a_tmp.obs[["n_genes_by_counts", "total_counts", "pct_counts_mt"]].copy()
        obs_tmp["Paper_ID"] = study
        all_qc_dfs.append(obs_tmp)
        a_tmp.file.close()
        continue

    study_meta  = meta_clean[meta_clean["Paper_ID"] == study]
    adata_list  = []

    for _, row in study_meta.iterrows():
        id_val = row["ID"]
        pid    = row["patient_id"]
        sid    = row["sample_id"]
        loc    = row["countmatrixlocation"]

        # Build path to Cell Ranger h5 output
        if study == "P02_2022_NatCommun_Wong":
            matrix_path = (
                f"/mctp/share/users/yupingz/cellranger/cellranger_7.1/counts/{loc}/"
                f"{pid}/outs/per_sample_outs/{pid}/count/sample_filtered_feature_bc_matrix.h5"
            )
        elif study == "P01_2021_CellRep_Tuong":
            matrix_path = (
                f"/mctp/share/users/hanbyul/cellranger/cellranger_7.1/counts/{loc}/"
                f"{sid}/outs/filtered_feature_bc_matrix.h5"
            )
        else:
            matrix_path = (
                f"/mctp/share/users/yupingz/cellranger/cellranger_7.1/counts/{loc}/"
                f"{sid}/outs/filtered_feature_bc_matrix.h5"
            )

        if not os.path.exists(matrix_path):
            print(f"  [missing] {matrix_path}")
            continue

        a = sc.read_10x_h5(matrix_path, gex_only=True)
        a.var_names_make_unique()

        # Wong 2022: subset to tumor or benign barcodes from HTO demultiplexing
        if study == "P02_2022_NatCommun_Wong":
            hto_label = "B" if row["Group"] == "B" else "T"
            wb = wb_by_patient[pid]
            valid_bc = wb.index[wb["final.call"] == hto_label]
            keep_bc  = a.obs_names.intersection(valid_bc)
            a = a[keep_bc].copy()

        # Prefix barcodes with sample ID (matches R convention)
        a.obs_names = [f"{id_val}_{bc}" for bc in a.obs_names]
        print(f"  {id_val}: {a.n_obs:,} cells loaded")
        adata_list.append(a)

    if not adata_list:
        print(f"  [warning] No samples loaded for {study}")
        continue

    # Merge samples for this study
    merged = sc.concat(adata_list, join="outer")
    merged.obs["Paper_ID"] = study

    # Compute QC metrics
    merged.var["mt"] = merged.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(merged, qc_vars=["mt"], inplace=True)

    print(f"  {study}: {merged.n_obs:,} cells, {merged.n_vars:,} genes")
    merged.write_h5ad(out_path)

    obs_qc = merged.obs[["n_genes_by_counts", "total_counts", "pct_counts_mt"]].copy()
    obs_qc["Paper_ID"] = study
    all_qc_dfs.append(obs_qc)

# Save pre-filter QC metrics (used by 05_visualization.R for Supplementary Figure 1a)
qc_obs = pd.concat(all_qc_dfs, axis=0)
qc_csv = os.path.join(PY_FINAL, "QC_metrics_raw_prefilter.csv")
qc_obs.to_csv(qc_csv)
print(f"\nSaved pre-filter QC metrics: {qc_csv}")
print(f"Total cells before any filtering: {len(qc_obs):,}")

# Concatenate all per-study h5ad into one merged object
print("\nConcatenating all studies...")
all_adata = []
for study in study_ids:
    path = os.path.join(DATA_DIR, f"{study}_combined.h5ad")
    if os.path.exists(path):
        all_adata.append(sc.read_h5ad(path))
    else:
        print(f"  [missing] {path}")

combined = sc.concat(all_adata, keys=study_ids, label="Paper_ID", join="outer", index_unique=None)
combined_path = os.path.join(DATA_DIR, "all_studies_combined.h5ad")
combined.write_h5ad(combined_path)
print(f"Saved: {combined_path}  ({combined.n_obs:,} cells, {combined.n_vars:,} genes)")

# =============================================================================
# Part B — Filter by R cell type annotation
# metacomb_17data_final.csv is the cell annotation table exported by
# 01_QC_and_annotation.R after SingleR + manual curation.
# Only cells present in this table are kept; all obs columns are replaced
# with the R annotations.
# =============================================================================
print("\n" + "="*60)
print("Part B: Filter by R annotation → merged_a17.h5ad")
print("="*60)

cellanno_path = os.path.join(BASE_DIR, "metacomb_17data_final.csv")
cellanno = pd.read_csv(cellanno_path, index_col=0, low_memory=False)
print(f"R annotation table: {len(cellanno):,} cells")

adata = sc.read_h5ad(combined_path)
print(f"Python merged: {adata.n_obs:,} cells")

# Keep only cells that exist in both the Python object and the R annotation
keep = cellanno.index.intersection(adata.obs_names)
print(f"Overlap (in both): {len(keep):,} cells")

adata = adata[keep].copy()
adata.obs = cellanno.loc[adata.obs_names].copy()   # replace obs with R annotations

merged_a17_path = os.path.join(PY_FINAL, "merged_a17.h5ad")
adata.write_h5ad(merged_a17_path)
print(f"Saved: {merged_a17_path}  ({adata.n_obs:,} cells)")

# =============================================================================
# Part C — Create integration inputs
# =============================================================================
print("\n" + "="*60)
print("Part C: Create integration inputs")
print("="*60)

# ── C1: All-cells integration input ──────────────────────────────────────────
# Remove low-quality / unannotated cells and mitochondrial genes.
# The column "final_hit1_group_broad" contains the broad cell type label from R.
adata = sc.read_h5ad(merged_a17_path)

mask_lq = adata.obs["final_hit1_group_broad"].isin(EXCLUDE_LABELS)
adata   = adata[~mask_lq].copy()
adata   = adata[:, ~adata.var_names.str.startswith("MT-")].copy()

# Store raw counts as a layer; log-normalize expression into X
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.layers["lognorm"] = adata.X.copy()

merged_hit1_path = os.path.join(PY_FINAL, "merged_a17_hit1.h5ad")
adata.write_h5ad(merged_hit1_path)
print(f"Saved all-cells input: {merged_hit1_path}  ({adata.n_obs:,} cells)")
print(f"  → used by 04_integration_all_cells.py")

# ── C2: Luminal cell subset for Seurat R batch correction ────────────────────
# Subset to Epithelial Cells only; further filter inside 03_batch_correction_luminal.R
# (No_Group luminal cells are removed there).
# Export a "stripped" h5ad (counts + lognorm only, no scale.data) because
# SeuratDisk / SaveH5Seurat fails with scale.data on large objects.
# Also export full obs metadata as a CSV so R can re-attach all columns.
lum_adata = adata[adata.obs["final_hit1_group_broad"] == "Epithelial Cells"].copy()
lum_adata = lum_adata[~lum_adata.obs["etsgroup"].isin({"No_Group", "", "nan"})].copy()
print(f"Luminal subset (No_Group removed): {lum_adata.n_obs:,} cells")

# Save full obs metadata as CSV for R to re-attach after h5ad → Seurat conversion
obs_csv_path = os.path.join(SCIB_DIR, "obs_metadata_lum_v1_noNoGroup.csv")
lum_adata.obs.to_csv(obs_csv_path)
print(f"Saved obs metadata CSV: {obs_csv_path}")

# Strip to bare h5ad (counts + lognorm layers only; no obsm, no varm, no uns)
lum_stripped = sc.AnnData(
    X    = lum_adata.layers["counts"].copy(),
    obs  = lum_adata.obs[[]].copy(),   # empty obs (metadata re-attached in R via CSV)
    var  = lum_adata.var[[]].copy(),
)
lum_stripped.layers["counts"]  = lum_adata.layers["counts"].copy()
lum_stripped.layers["lognorm"] = lum_adata.layers["lognorm"].copy()
lum_stripped.obs_names = lum_adata.obs_names
lum_stripped.var_names = lum_adata.var_names

stripped_path = os.path.join(SCIB_DIR, "merged_a17_STRIPPED_lum_v1_noNoGroup.h5ad")
lum_stripped.write_h5ad(stripped_path)
print(f"Saved stripped luminal h5ad: {stripped_path}")
print(f"  → used by 03_batch_correction_luminal.R (CCA / RPCA / Harmony)")

print("\nPart C complete.")
print("Next steps:")
print("  1. Run 03_batch_correction_luminal.R  (CCA + RPCA + Harmony in Seurat v5)")
print("  2. Run 04_integration_all_cells.py    (Harmony + Scanorama + BBKNN + scVI + scANVI)")
print("  3. Run 05_integration_luminal_cells.py (imports R embeddings + Scanorama + BBKNN + scVI + scANVI)")
