# A harmonized single-cell RNA-seq atlas of human localized and metastatic prostate cancers and benign tissues: code

This repository holds the code of the main processing steps that made the atlas files deposited at Zenodo. Without the steps and inputs listed under Not included, the scripts do not run end to end.

## Article, data and citation

- Article: Data Descriptor submitted to Scientific Data.
- Data: Zenodo, https://doi.org/10.5281/zenodo.19500694.
- Code: release v1.0 of this repository, archived at Zenodo.
- Repository: https://github.com/hanbyulcho/prostate-scrna-atlas

If you use this code or the data, cite the article and the Zenodo deposits of the data and code.

## Terms

- Barcode-swap studies: the four studies (47 samples) in which swappedDrops removed barcode-swapped molecules.
- Re-call: the second cell calling of those studies by emptyDrops, on the counts left after that removal.
- STUDY names of the re-call scripts: P01_Tuong, P06_Heidegger, P07_Chan and P09_Cheng stand for the study_id values P01_2021_CellRep_Tuong, P06_2022_MolCancer_Heidegger, P07_2022_Science_Chan and P09_2022_EurUrol_Cheng of sample_metadata.csv.
- First atlas version: the atlas built before the re-call. It is not deposited, and several scripts read its cell metadata.
- Atlas cell types: the 18 cell types of the atlas, stored as cell_type (author_cell_type in the CELLxGENE-schema versions).
- Source labels: the labels the original papers gave.
- Curator: the author who decided the labels of flagged clusters, of removed cell groups and of luminal units.
- Restored cells: cells of a group (setB_unit in the code) that the first atlas version removed and the curator restored, with one label per group.
- Luminal unit: the kept luminal cells of one cluster in one sample.
- Luminal subgroup label: the putative label that the labeling rule or the curator gave each luminal unit, stored as luminal_subgroup.

## Scripts

The table lists the scripts in `pipeline/` that are run directly, in run order. Each script comes after the scripts whose outputs it reads. Every script in `pipeline/` except `verbatim_functions.R` opens with a header note. In `umap_sigma_fix.py` the note states what the patch does and why, and in the other scripts it covers their inputs, key settings and Methods subsections.

| Script | Task | Environment |
|---|---|---|
| `02_qc/swap/check_barcode_sharing.py` | Tests excess barcode sharing between GEM wells | Python 3 with h5py and scipy (versions not recorded) |
| `02_qc/swap/04_swappeddrops_study.R` | Removes barcode-swapped molecules | `r433` |
| `03_annotation/per_sample/01_qc_doublets.R` | Applies cell QC and calls doublets in one sample, for L071, L111 and L126 in the run | `r433` |
| `03_annotation/per_sample/02_cluster_annotate.R` | Clusters one sample and labels the clusters with SingleR, for L071, L111 and L126 in the run (loads `verbatim_functions.R` and `helpers.R`) | `r433` |
| `02_qc/recall/03_recall_qc_classify.R` | Applies the emptyDrops re-calls and cell QC | `r433` |
| `02_qc/recall/04_doublets_cluster_annotate.R` | Calls doublets, clusters and labels the re-called cells of each sample of a barcode-swap study | `r433` |
| `02_qc/recall/04b_restored_cells.R` | Same, for those samples with restored cells | `r433` |
| `02_qc/soupx/01_soupx_unit.R` | Estimates the ambient-RNA contamination fraction (SoupX), without correcting counts | `r433` |
| `08_copy_number/arms/01_pseudobulk_v4.R` | Sums raw counts into pseudobulks of luminal units and of reference cell groups | `r433` |
| `08_copy_number/arms/02_cnv_v4.R` | Infers arm-level copy number | `r433` |
| `08_copy_number/arms/NG_L5_v2.R` | Same, for luminal cells outside the first atlas version | `r433` |
| `07_luminal_labels/R_rules_v9.R` | Applies the luminal labeling rule | `r433` |
| `05_integration/allcells/03_run_python_integration.py` | Integrates all cells and runs scib-metrics (loads `umap_sigma_fix.py`) | `python311` |
| `05_integration/luminal/02_r_seurat_luminal.R` | Runs Seurat CCA, RPCA and Harmony | `r440_seurat` |
| `05_integration/luminal/03_run_python_integration.py` | Integrates luminal cells and runs scib-metrics (loads `umap_sigma_fix.py`) | `python311` |
| `07_luminal_labels/label_transfer/S1_soft_knn_primary_v1.py` | Transfers labels by scArches weighted kNN (loads `lu_common_v1.py`) | `python311` |
| `06_benchmark/R1_rescore_scib.py` | Scores the stored luminal embeddings on the luminal subgroup labels, without re-running integration (loads `umap_sigma_fix.py` from the script folder of the stored integration run). The means of its per-seed scores are in the deposit's `scib_luminal_raw.csv` | `python311` |
| `06_benchmark/04_scib_tables_figures.py` | Makes the score tables and diagnostic plots | `python311` |
| `13_cellxgene/convert_to_cellxgene.py` | Converts the all-cell and luminal objects to CELLxGENE schema 7.1.0 (loads `cxg_mappings.py`) | `cxg` |

`verbatim_functions.R` holds sc.proc, get.psk and process_singleR, copied unchanged from the code of the first atlas version. Sheet1 of `celltype_colorcode_update_17papers.xlsx` maps each SingleR reference label (column detail) to a reduced label (detail_reduce) and a broader type (group_detailed). `helpers.R` turns the group_detailed values into 17 of the 18 atlas cell types. The curator's cluster decisions (not included) add Ionocyte, the 18th atlas cell type. `umap_sigma_fix.py` corrects the sigma = inf that umap-learn 0.5.9.post2 with numba 0.62.1 gives, so the UMAPs, BBKNN graphs and kBET scores stay correct.

For release, we edited only comments, docstrings, messages, record texts, output notes, xlsx file properties, some file names, some paths and one check's patterns. The record texts are the values of seven constants in `cxg_mappings.py` that name the source of a decided value. Output notes are the text values that the scripts write into output tables and run settings, such as the reason column and the seed_scheme notes. In them, we replaced decision codes, internal file names and the curator's name with plain words. The status text of swap-step removals changed, but no cell is kept, removed or labeled differently. In the scripts, each path to a folder or file of the authors' system begins with a placeholder. In the five R scripts that load a 00_env.R file (not included) and in helpers.R, other paths begin with folders that this file sets. One more path edit makes both `03_run_python_integration.py` scripts load `umap_sigma_fix.py` from `05_integration/`. In the public-text check of `cxg_mappings.py`, the path pattern /home/ replaced three personal-name patterns.

The path placeholders stand for these files and folders:

- `<input: ...>`, `<input folder: ...>`, `<work folder: ...>` and `<output folder: ...>`: the file or folder that the placeholder's text after its colon describes
- `<PROJECT_DIR>`: the project folder that these scripts read from and write to
- `<BASE_DIR>`: the folder whose 17datasets/ subfolder holds files of the first atlas version
- `<R_LIBRARY>`: the main package library of the `r433` environment
- `<R_HOME>`: the R 4.3.3 installation folder
- `<PY311_ENV>`: the folder of the `python311` virtual environment (human_sc in the restore commands below)
- `<CXG_ENV>`: the folder of the `cxg-schema` conda environment that `env/cxg/environment.yml` records
- `<LEIDEN_ENV>`: the folder of the leiden-fixed conda environment, not recorded in `env/` (leidenalg 0.10.2, python-igraph 0.11.8), whose lib/ folder goes on LD_LIBRARY_PATH for 04l_lumrerun_child.R
- `<SCRATCH_DIR>`: a scratch folder for temporary and work files
- `<CELLRANGER_RERUN_DIR>`: a folder of Cell Ranger re-runs, which `check_barcode_sharing.py` reads for L111 and L126

Scripts that run in `r433` search `<PROJECT_DIR>/redo_env_v1/Rlib_R433_overlay` before `<R_LIBRARY>`, so they use the `env/r433/renv.lock` versions of its 42 packages.

`04_doublets_cluster_annotate.R` and `04b_restored_cells.R` stop unless `verbatim_functions.R` has the md5 of the run copy, which had a two-line header. To rerun them with the released file, replace 35609a43c33765adf8123294db7914da with 568d299a151c8e1b1e0a5ed9ed189c41 in both scripts.

## Not included

The steps and inputs below are not in this repository. Some items state why.

- the Cell Ranger runs and outputs, the sample manifests and the molecule_info.h5 files
- the five SingleR reference objects (prostateref, metref, p1, m2 and m2_group in the code). The prostateref object is the article's annotation reference, and these scripts label each cluster with its best-matching prostateref label.
- the objects and cell metadata of the first atlas version
- the per-sample processing of the first atlas version: cell QC, the scDblFinder run recorded as 'run 1' in the deposit, clustering and SingleR labeling
- 01_prepare_inputs.py and its outputs
- 00_env.R, libraries.R and basic.info.rds
- step 1 of the re-call, which made the swap-cleaned counts of all barcodes of each study
- the emptyDrops check at the DropletUtils defaults and its re-run for some samples (step 2), whose calls the re-call uses
- three steps for some samples only, left out to keep the repository small:
  - the per-GEM-well processing of L138 and L139
  - the epithelial re-clustering in 14 samples, including 04l_lumrerun_child.R
  - the L132 rule that keeps a barcode only if it passes cell QC on both of its Cell Ranger runs
- the evidence step of the luminal labeling rule (unit scores, ETS tests against the SoupX ambient-RNA expectation and support vector machine classes), which R_rules_v9.R reads
- the averaging of seeds 0 and 1 for the score tables
- the relabeling of 141 cells of sample L082 from T cell to Luminal after both integrations
- the assembly of qc_summary_sample_level.csv
- metadata curation, gene-set selection, the comparison with the source labels and the other CELLxGENE steps, such as writing gene_id_remapping.csv and cellxgene_field_mapping.csv
- the code of the article figures, left out because the repository holds processing code only
- applying the curator's decisions and writing the deposited objects and tables from them. These scripts read the curator's decision records, which are not deposited. master_cell_table.csv.gz, all_cells_metadata.csv and luminal_cells_metadata.csv record each cell's result.
- any other step or input that the scripts name and this repository does not hold

## Software environments

`env/` records four environments, and `env/versions_summary.tsv` lists their key versions. Cell Ranger 7.1.0 (GRCh38-2020-A) made the count matrices, after bamtofastq (1.4.1 where recorded) for BAM input. We did not test the restore commands below.

```bash
# env/python311 (Python 3.11.1). torch 2.5.1+cpu is from the PyTorch wheel index.
python3.11 -m venv human_sc && . human_sc/bin/activate
python -m pip install --no-deps --index-url https://download.pytorch.org/whl/cpu torch==2.5.1+cpu
python -m pip install --no-deps -r env/python311/requirements.txt
# env/cxg (Python 3.11.17, cellxgene-schema 7.1.1)
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

Table 2 of the article lists the files of the Zenodo deposit. The code below opens the main files with packages that `env/python311` and `env/r433` pin.

In `all_cells_integrated.h5ad` and `luminal_cells_integrated.h5ad`, `X` holds the raw counts. Their `obsm` holds, for seeds 0 and 1, each benchmarked method's embedding and UMAP coordinates (for BBKNN, UMAP coordinates only), such as `X_scanorama_seed0` and `X_umap_scanorama_seed0`. To redraw the article's UMAPs, use `X_umap_scanvi_seed0` for all cells (Figures 2a and 5c) and `X_umap_scanorama_seed0` for luminal cells (Figures 5d,e and 6a–c).

In the CELLxGENE-schema versions (`all_cells_cellxgene_v1.h5ad` and `luminal_cells_cellxgene_v1.h5ad`), `X` holds log-normalized expression and `raw.X` the raw counts. Their atlas cell types are in `author_cell_type`, and `cell_type_ontology_term_id` holds the Cell Ontology term. The integrated objects index genes by symbol, and the CELLxGENE-schema versions index them by GENCODE v48 Ensembl identifier.

Python, with anndata 0.12.4 and pandas 2.3.3:

```python
import anndata as ad
import pandas as pd

lum = ad.read_h5ad("luminal_cells_integrated.h5ad")
print(list(lum.obsm))  # embedding and UMAP keys
allc = ad.read_h5ad("all_cells_integrated.h5ad", backed="r")  # X stays on disk
cxg = ad.read_h5ad("luminal_cells_cellxgene_v1.h5ad")
counts = cxg.raw.X

cells = pd.read_csv("luminal_cells_metadata.csv")
samples = pd.read_csv("sample_metadata.csv")
cells = cells.merge(samples[["library_id", "castration_status", "metastatic_status"]],
                    on="library_id", how="left", validate="many_to_one")  # sample fields not in the cell table
qc = pd.read_csv("qc_metrics_prefilter.csv.gz", low_memory=False)  # column types from all rows
```

In R, `reader = "R"` selects the R reader of zellkonverter, which needs no Python. That reader does not read `raw.X`, so in R take the raw counts from the integrated objects.

R, with SingleCellExperiment 1.24.0, zellkonverter 1.17.0, Seurat 5.1.0, data.table 1.15.4 and readxl 1.4.3:

```r
library(SingleCellExperiment)
library(zellkonverter)
library(data.table)

lum <- readH5AD("luminal_cells_integrated.h5ad", reader = "R")  # a SingleCellExperiment
assayNames(lum)  # "X": raw counts, genes x cells
reducedDimNames(lum)  # embedding and UMAP keys
allc <- readH5AD("all_cells_integrated.h5ad", reader = "R", use_hdf5 = TRUE)  # counts stay on disk
seu <- Seurat::as.Seurat(lum, counts = "X", data = NULL)  # assay "originalexp", one key-renaming warning per reduction

cells <- fread("luminal_cells_metadata.csv")
samples <- fread("sample_metadata.csv")
cells <- merge(cells, samples[, c("library_id", "castration_status", "metastatic_status")],
               by = "library_id", all.x = TRUE)  # sample fields not in the cell table
qc <- fread("qc_metrics_prefilter.csv.gz")  # reading .gz needs R.utils, which env/r433 pins
per_study <- readxl::read_excel("original_annotation_comparison.xlsx", sheet = "Per_study")
```

The deposit's README.md gives the format, join keys and conventions of its tables. Table 2 names the tables whose columns `metadata_dictionary.csv` describes. A spreadsheet program also opens the CSV files, after decompressing a .csv.gz file (for example with `gunzip -k`). Excel shows at most 1,048,576 rows, so read `qc_metrics_prefilter.csv.gz` (1,182,729 rows) with pandas or R.

## License

MIT (`LICENSE`). `lu_common_v1.py` also keeps the scArches BSD 3-Clause license notice, because one of its checks adapts the kernel-weight code of scArches 0.6.1 (scarches/utils/knn.py).

## Contact

Open an issue here, or write to the article's corresponding author.
