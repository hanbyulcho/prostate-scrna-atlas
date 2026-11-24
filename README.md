# A harmonized single-cell RNA-seq atlas of human localized and metastatic prostate cancers and benign tissues: code

This repository holds the code of the main processing steps that made the atlas files deposited at Zenodo. The scripts do not run end to end without the steps and inputs listed under [Not included](#not-included).

## Article, data and citation

- Article: Data Descriptor submitted to Scientific Data.
- Data: Zenodo, https://doi.org/10.5281/zenodo.19500694.
- Code: release v1.0 of this repository, archived at Zenodo.
- Repository: https://github.com/hanbyulcho/prostate-scrna-atlas

If you use this code or the data, cite the article and both Zenodo deposits.

## Scripts

The table lists the scripts in `pipeline/` that are run directly, in run order. Every script in `pipeline/` except `verbatim_functions.R` opens with a header note. Environment names refer to folders in `env/`.

The re-call is a second cell calling by emptyDrops in four studies (47 samples), on the counts left after swappedDrops removed barcode-swapped molecules. The first atlas version, built before the re-call, is not deposited. A luminal unit is the luminal cells of one cluster in one sample. The curator is the author who decided the labels of flagged clusters, of removed cell groups and of luminal units.

| Script | Task | Environment |
|---|---|---|
| `02_qc/swap/check_barcode_sharing.py` | Tests excess barcode sharing between GEM wells | Python 3 with h5py and scipy (versions not recorded) |
| `02_qc/swap/04_swappeddrops_study.R` | Removes barcode-swapped molecules | `r433` |
| `03_annotation/per_sample/01_qc_doublets.R` | Applies cell QC and calls doublets in one sample (run for L071, L111 and L126) | `r433` |
| `03_annotation/per_sample/02_cluster_annotate.R` | Clusters one sample and labels the clusters with SingleR, for L071, L111 and L126 (loads `verbatim_functions.R` and `helpers.R`) | `r433` |
| `02_qc/recall/03_recall_qc_classify.R` | Applies the re-calls and cell QC | `r433` |
| `02_qc/recall/04_doublets_cluster_annotate.R` | Calls doublets, then clusters and labels the re-called cells of each sample | `r433` |
| `02_qc/recall/04b_restored_cells.R` | Same, for samples with cells that the first atlas version removed and the curator restored | `r433` |
| `02_qc/soupx/01_soupx_unit.R` | Estimates the ambient-RNA fraction (SoupX) without correcting counts | `r433` |
| `08_copy_number/arms/01_pseudobulk_v4.R` | Sums raw counts into pseudobulks of luminal units and of reference cell groups | `r433` |
| `08_copy_number/arms/02_cnv_v4.R` | Infers arm-level copy number | `r433` |
| `08_copy_number/arms/NG_L5_v2.R` | Same, for luminal cells outside the first atlas version | `r433` |
| `07_luminal_labels/R_rules_v9.R` | Applies the luminal labeling rule | `r433` |
| `05_integration/allcells/03_run_python_integration.py` | Integrates all cells and runs scib-metrics (loads `umap_sigma_fix.py`) | `python311` |
| `05_integration/luminal/02_r_seurat_luminal.R` | Runs Seurat CCA, RPCA and Harmony | `r440_seurat` |
| `05_integration/luminal/03_run_python_integration.py` | Integrates luminal cells and runs scib-metrics (loads `umap_sigma_fix.py`) | `python311` |
| `07_luminal_labels/label_transfer/S1_soft_knn_primary_v1.py` | Transfers labels by scArches weighted kNN (loads `lu_common_v1.py`) | `python311` |
| `06_benchmark/R1_rescore_scib.py` | Scores the stored luminal embeddings on the deposited `luminal_subgroup` labels (loads `umap_sigma_fix.py`). The means of its seed 0 and seed 1 scores are in the data deposit's `scib_luminal_raw.csv` | `python311` |
| `06_benchmark/04_scib_tables_figures.py` | Makes the score tables and diagnostic plots | `python311` |
| `13_cellxgene/convert_to_cellxgene.py` | Converts the all-cell and luminal objects to CELLxGENE schema 7.1.0 (loads `cxg_mappings.py`) | `cxg` |

`verbatim_functions.R` holds three functions copied unchanged from the first atlas version. `helpers.R` and Sheet1 of `celltype_colorcode_update_17papers.xlsx` turn SingleR labels into atlas cell types. `umap_sigma_fix.py` corrects a sigma = inf from umap-learn 0.5.9.post2 with numba 0.62.1. Uncorrected, this value would distort the UMAPs, BBKNN graphs and kBET scores.

For release, we edited only text, such as comments, messages, file names, paths and text written into outputs. These edits change no cell's label and do not change which cells are kept or removed.

Paths to the authors' system begin with a placeholder, such as `<PROJECT_DIR>`, or with a folder set in a `00_env.R` file (not included). `<LEIDEN_ENV>` is the folder of a conda environment with leidenalg 0.10.2 and python-igraph 0.11.8, which `env/` does not record.

`04_doublets_cluster_annotate.R` and `04b_restored_cells.R` check the md5 of `verbatim_functions.R`. To run them with the released file, replace 35609a43c33765adf8123294db7914da with 568d299a151c8e1b1e0a5ed9ed189c41 in both.

## Not included

Steps and inputs not in this repository include the following.

- Cell Ranger runs and outputs, the five SingleR reference objects and the first atlas version (objects, cell metadata and per-sample processing)
- `00_env.R`, `libraries.R`, `basic.info.rds` and `01_prepare_inputs.py`
- re-call step 1, the emptyDrops check at DropletUtils defaults and step 2 (its rerun for some samples)
- other steps run for some samples only, such as the epithelial re-clustering in 14 samples
- the evidence step of the luminal labeling rule, which `R_rules_v9.R` reads
- applying the curator's decisions and writing the deposited integrated objects and tables, such as the seed-averaged score tables
- metadata curation, gene-set selection, the comparison with the original papers' labels and the other CELLxGENE steps
- the code of the article figures

## Software environments

`env/` records four environments, and `env/versions_summary.tsv` lists their key versions. Cell Ranger 7.1.0 (GRCh38-2020-A) made the count matrices, after bamtofastq (1.4.1 where recorded) for BAM input. We did not test the restore commands below.

```bash
# env/python311 (Python 3.11.1). torch 2.5.1+cpu is from the PyTorch wheel index.
python3.11 -m venv human_sc && . human_sc/bin/activate
python -m pip install --no-deps --index-url https://download.pytorch.org/whl/cpu torch==2.5.1+cpu
python -m pip install --no-deps -r env/python311/requirements.txt
# env/cxg (Python 3.11.17, cellxgene-schema package 7.1.1)
conda env create --file env/cxg/environment.yml
```

```r
install.packages("renv", repos = "https://cloud.r-project.org")
# in env/r433/ with R 4.3.3
renv::restore(lockfile = "renv.lock", exclude = "ifnb.SeuratData")
# in env/r440_seurat/ with R 4.4.0
renv::restore(lockfile = "renv.lock", exclude = c("bluster", "BayesPrism", "panc8.SeuratData"))
```

renv cannot restore the excluded packages as recorded, and these scripts do not load them. The `env/r433` lockfile does not list Matrix, which comes with R 4.3.3.

## Loading the deposited files

The data deposit's README.md describes each file and how to read its tables. In the two `*_integrated.h5ad` files, `X` holds the raw counts and genes are indexed by symbol. Their `obsm` holds each benchmarked method's embedding and UMAP coordinates for seeds 0 and 1 (for BBKNN, UMAP coordinates only). To redraw the article's UMAPs, use `X_umap_scanvi_seed0` for all cells and `X_umap_scanorama_seed0` for luminal cells.

In the two `*_cellxgene_v1.h5ad` files, `X` holds log-normalized expression and `raw.X` the raw counts. Their genes are indexed by GENCODE v48 Ensembl identifier, and the atlas cell types are in `author_cell_type`.

Python, with anndata 0.12.4:

```python
import anndata as ad

lum = ad.read_h5ad("luminal_cells_integrated.h5ad")
print(list(lum.obsm))  # embedding and UMAP keys
allc = ad.read_h5ad("all_cells_integrated.h5ad", backed="r")  # X stays on disk
```

The zellkonverter option `reader = "R"` needs no Python but does not read `raw.X`, so in R take raw counts from the integrated objects.

R, with SingleCellExperiment 1.24.0 and zellkonverter 1.17.0:

```r
library(SingleCellExperiment)
library(zellkonverter)

lum <- readH5AD("luminal_cells_integrated.h5ad", reader = "R")  # assay "X": raw counts
reducedDimNames(lum)  # embedding and UMAP keys
allc <- readH5AD("all_cells_integrated.h5ad", reader = "R", use_hdf5 = TRUE)  # counts stay on disk
```

## License

MIT (`LICENSE`). `lu_common_v1.py` also keeps the scArches BSD 3-Clause license notice for code adapted from scArches 0.6.1.

## Contact

Open an issue in this repository or contact Hanbyul Cho (@hanbyulcho).
