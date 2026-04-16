# A Harmonized Single-Cell RNA-seq Atlas of Human Prostate Cancer

**Cho et al.** — *Scientific Data* (in submission)

This repository contains the complete analysis code for integrating 17 published prostate single-cell RNA-seq datasets into a harmonized atlas spanning benign tissue, primary prostate cancer (PCa), and metastatic castration-resistant prostate cancer (mCRPC).

---

## Atlas Overview

| | |
|---|---|
| **Datasets** | 17 studies |
| **Samples** | 163 |
| **Donors** | 106 |
| **Final cells** | 754,195 |
| **Integration method** | scANVI |
| **Disease groups** | Benign (N/B/BPH), Primary PCa, mCRPC |

---

## Pipeline

The analysis runs in 7 sequential steps across R and Python:

```
01_QC_and_annotation.R          (R)
        ↓
02_data_merge_and_filter.py     (Python)
        ↓
03_batch_correction_luminal.R   (R)
        ↓
04_integration_all_cells.py     (Python)
05_integration_luminal_cells.py (Python)
        ↓
06_h5ad_to_seurat.R             (R)
        ↓
07_visualization.R              (R)
```

---

## Script Descriptions

### `01_QC_and_annotation.R`
**Language:** R | **Key packages:** Seurat v5, scDblFinder, SingleR

- Loads raw 10x Genomics count matrices (Cell Ranger v7.1 output) for each of the 17 studies
- Applies per-sample QC filters: ≥250 genes, ≥500 UMIs, <20% mitochondrial reads
- Detects doublets with **scDblFinder**
- Generates pseudobulk profiles per cluster and annotates cell types using **SingleR** with a prostate-specific reference (CellAnno + Henry et al.)
- Exports `all_cells_prefilter_annotation.csv` — the master cell annotation table used by all downstream steps

### `02_data_merge_and_filter.py`
**Language:** Python | **Key packages:** scanpy, anndata

- Loads raw 10x `.h5` files per sample into AnnData objects (includes special HTO demultiplexing for Wong 2022 / P02)
- Computes QC metrics (nCount, nFeature, pct_mt) and saves `QC_metrics_prefilter.csv`
- Concatenates all samples → `all_studies_combined.h5ad`
- Merges R annotation from `all_cells_prefilter_annotation.csv`, filters to annotated cells → `merged_a17.h5ad`
- Removes low-quality / unknown cell types and MT genes → **`all_cells_raw.h5ad`** (main all-cells input)
- Subsets to luminal/epithelial cells → stripped h5ad + obs CSV for R batch correction

### `03_batch_correction_luminal.R`
**Language:** R | **Key packages:** Seurat v5, SeuratDisk

Runs Seurat v5 `IntegrateLayers()` for three methods on luminal cells:

| Method | Seurat call | Output reduction |
|--------|------------|-----------------|
| **CCA** | `CCAIntegration` | `integrated.cca` / `umap.cca` |
| **RPCA** | `RPCAIntegration` | `integrated.rpca` / `umap.rpca` |
| **Harmony** | `HarmonyIntegration` | `harmony` / `umap.harmony` |

Uses `DietSeurat` + `SaveH5Seurat` + `Convert` to export **`luminal_cells_raw.h5ad`** — a single h5ad carrying all three sets of embeddings alongside raw/lognorm expression.

### `04_integration_all_cells.py`
**Language:** Python | **Key packages:** scanpy, scvi-tools, scib-metrics

- Input: `all_cells_raw.h5ad` (754,195 cells)
- Parameters: **2,000 HVGs, 50 PCs**, `batch_key = Paper_ID_short`, `label_key = final_hit1_group_detailed`
- Runs and benchmarks 6 integration methods:

| Method | Type | scib Total |
|--------|------|-----------|
| Unintegrated (PCA) | Baseline | 0.48 |
| Harmony | Embedding | 0.69 |
| Scanorama | Embedding | 0.40 |
| BBKNN | Graph | 0.33 |
| scVI | Deep learning | 0.41 |
| **scANVI** | Deep learning | **0.74** ✓ |

- Exports scib benchmark tables → **Supplementary Tables 4 & 5**
- Saves `all_cells_integrated_stripped.h5ad` + `all_cells_metadata.csv`

### `05_integration_luminal_cells.py`
**Language:** Python | **Key packages:** scanpy, scvi-tools, scib-metrics

- Input: `luminal_cells_raw.h5ad` (from step 03, carrying CCA/RPCA/Harmony embeddings)
- Parameters: **3,000 HVGs, 30 PCs**, `batch_key = Paper_ID_short`, `label_key = etsgroup`
- Step 1: Recovers lognorm from Seurat's raw slot; back-computes pseudo-counts for scVI
- Step 2: Renames Seurat embedding keys (e.g. `X_integrated.cca` → `X_cca`)
- Benchmarks all **8 methods** (4 from R + 4 from Python):

| Method | Source | scib Total |
|--------|--------|-----------|
| Unintegrated | Python | 0.46 |
| CCA | R (Seurat) | 0.36 |
| RPCA | R (Seurat) | 0.41 |
| Harmony | R (Seurat) | 0.51 |
| Scanorama | Python | 0.46 |
| BBKNN | Python | 0.26 |
| scVI | Python | 0.54 |
| **scANVI** | Python | **0.78** ✓ |

- Luminal ETS subgroups: Benign, BPH, ERG, ETV, nonETS, AR, NEPC, DN (AR−/NE−)
- Exports scib benchmark tables → **Supplementary Tables 6 & 7**

### `06_h5ad_to_seurat.R`
**Language:** R | **Key packages:** Seurat v5, SeuratDisk

Converts the Python scib output h5ads back to Seurat RDS for R-based visualization:
1. `Convert()` — h5ad → h5seurat (SeuratDisk intermediate format)
2. `LoadH5Seurat()` — load into Seurat object
3. `AddMetaData()` — re-attach full obs metadata from companion CSV
4. `UpdateSeuratObject()` + `saveRDS()`

Output: `all_cells_integrated.rds` and `luminal_cells_integrated.rds`

### `07_visualization.R`
**Language:** R | **Key packages:** Seurat v5, ggplot2, patchwork, ggrastr

---

## Dependencies

### R packages
```r
install.packages(c("Seurat", "ggplot2", "dplyr", "tidyr", "patchwork", "gridExtra"))
remotes::install_github("mojaveazure/seurat-disk")   # SeuratDisk
BiocManager::install(c("scDblFinder", "SingleR", "celldex"))
install.packages("ggrastr")
```

### Python packages
```bash
pip install scanpy anndata scvi-tools scib-metrics
```

---

## Data Availability

> **Note:** Data deposit links will be added upon paper acceptance.

Raw count matrices and processed AnnData/Seurat objects will be deposited at:
- **Zenodo**: *(link to be added upon acceptance)*

---

## Citation

Cho H, et al. A harmonized single-cell RNA-seq atlas of the human localized and metastatic prostate cancers and benign tissues. *Scientific Data* (2026, in submission).
