# Code release header: cell QC and doublet calls for one sample (one GEM well). In the run, it processed three samples: L071, L111 and L126.
# Input: the sample's count matrix (Cell Ranger filtered matrix or swap-cleaned counts), set in libraries.R (not included).
# Output: <LIB>_step01_qc_doublets.rds, which 02_cluster_annotate.R reads.
# Cell QC: at least 250 genes and 500 UMIs, and under 20% mitochondrial UMIs. A sample needs at least 100 passing cells.
# The 13 MT- genes are then removed. scDblFinder runs at its defaults after set.seed(SEED_DBL), with SEED_DBL = 0 (00_env.R, not included).
# scDblFinder is stochastic, and the released doublet classes are those in the doublet_class column of qc_metrics_prefilter.csv.gz. Methods: Quality control.
# 01_qc_doublets.R <LIB> -- cell QC, library rule, MT-gene removal, scDblFinder. One GEM well (capture) per library.
#   Read10X_h5 -> cell names "<LIB>_<barcode>" -> CreateSeuratObject, NormalizeData, percent.mt / ribo / hb (PercentageFeatureSet)
#   -> library rule of at least 100 QC-pass cells (a second rule, more than 100 barcodes with >= 250 genes, is only recorded)
#   -> QC-pass cells (metrics on all genes) -> the 13 MT genes removed
#   -> CreateSeuratObject, NormalizeData, FindVariableFeatures, ScaleData, as.SingleCellExperiment, set.seed(SEED_DBL), scDblFinder()
#   -> scDblFinder score and class saved -> singlets = class "singlet"
source("<PROJECT_DIR>/new_libraries_v1/scripts/00_env.R")
source(file.path(OUT, "scripts/libraries.R"))
suppressPackageStartupMessages({library(Seurat); library(Matrix); library(scDblFinder); library(SingleCellExperiment); library(dplyr)})
LIB <- commandArgs(TRUE)[1]; L <- get_lib(LIB)
env_report(); check_disk()
cat("library:", LIB, "| h5:", L$h5, "| count source:", L$count_source, "\n")
cat("h5 md5:", tools::md5sum(L$h5), " size:", file.size(L$h5), " mtime:", format(file.mtime(L$h5)), "\n")
dir.create(file.path(SCR, LIB), showWarnings = FALSE, recursive = TRUE)

basic.info <- readRDS(file.path(C, "basic.info.rds"))       # element 2 = the 13 MT- genes, plus $ribo and $hb_genes
mt <- basic.info[[2]]; ribo_genes <- basic.info$ribo; hb_genes <- basic.info$hb_genes
cat("MT genes (", length(mt), "):", mt, "\n"); stopifnot(length(mt) == 13, all(grepl("^MT-", mt)))

cts <- Read10X_h5(L$h5)
if (is.list(cts)) { cat("h5 modalities:", names(cts), "-> using 'Gene Expression'\n"); cts <- cts$`Gene Expression` }
cat("matrix:", nrow(cts), "genes x", ncol(cts), "barcodes; class", class(cts), "\n")
stopifnot(all(grepl("-1$", colnames(cts))), !any(duplicated(colnames(cts))))
cat("all 13 MT genes in the matrix:", all(mt %in% rownames(cts)), "; other rows starting 'MT-':", setdiff(grep("^MT-", rownames(cts), value = TRUE), mt), "\n")
stopifnot(all(mt %in% rownames(cts)))
cat("non-integer counts:", sum(cts@x != round(cts@x)), "; zero-UMI barcodes:", sum(Matrix::colSums(cts) == 0), "\n")
raw_bc <- colnames(cts)
colnames(cts) <- paste(LIB, raw_bc, sep = "_")

d <- CreateSeuratObject(counts = cts)
d <- NormalizeData(object = d, verbose = FALSE)
d[["percent.mt"]]   <- PercentageFeatureSet(d, features = mt)
d[["percent.ribo"]] <- PercentageFeatureSet(d, features = ribo_genes)
d[["percent.hb"]]   <- PercentageFeatureSet(d, features = hb_genes)

cells <- data.frame(cell = colnames(d), barcode = raw_bc, library = LIB, gem_well = L$gem_well,
                    nCount_RNA_allgenes = d$nCount_RNA, nFeature_RNA_allgenes = d$nFeature_RNA,
                    pct_mt = d$percent.mt, pct_ribo = d$percent.ribo, pct_hb = d$percent.hb, stringsAsFactors = FALSE)
cells$pass_nFeature250 <- cells$nFeature_RNA_allgenes >= 250
cells$pass_nCount500   <- cells$nCount_RNA_allgenes >= 500
cells$pass_mt20        <- cells$pct_mt < 20
cells$pass_qc <- cells$pass_nFeature250 & cells$pass_nCount500 & cells$pass_mt20      # cell QC (metrics on all genes)
n_ge250 <- sum(cells$pass_nFeature250); n_pass <- sum(cells$pass_qc)
lib_pass_original_rule <- n_ge250 > 100        # recorded only: more than 100 barcodes with >= 250 genes
lib_pass <- n_pass >= 100                      # library rule (applied): at least 100 QC-pass cells
cells$library_pass <- lib_pass
cat(sprintf("Cell Ranger-called barcodes %d | nFeature>=250 %d | nCount>=500 %d | MT<20 %d | pass all %d\n",
            nrow(cells), n_ge250, sum(cells$pass_nCount500), sum(cells$pass_mt20), n_pass))
cat("library rule (>=100 QC-pass cells):", lib_pass, "| rule of the first atlas version (>100 cells with >=250 genes):", lib_pass_original_rule, "\n")

out <- list(library = LIB, h5 = L$h5, h5_md5 = unname(tools::md5sum(L$h5)), cells = cells, n_ge250 = n_ge250, n_pass = n_pass,
            lib_pass = lib_pass, lib_pass_original_rule = lib_pass_original_rule, seed = SEED_DBL)
if (lib_pass) {
  cf <- cts[, cells$cell[cells$pass_qc], drop = FALSE]
  cf <- cf[!(rownames(cf) %in% c(mt)), , drop = FALSE]      # remove the 13 MT genes
  seu <- CreateSeuratObject(counts = cf)
  seu <- NormalizeData(object = seu)
  seu <- FindVariableFeatures(object = seu)
  seu <- ScaleData(object = seu)
  sce <- as.SingleCellExperiment(seu)
  cat("SCE assays:", assayNames(sce), "| dims", dim(sce), "\n")
  set.seed(SEED_DBL)                                         # fixed seed (SEED_DBL = 0, 00_env.R)
  sce <- scDblFinder(sce)                                    # default settings
  meta_scdblfinder <- sce@colData@listData %>% as.data.frame() %>% dplyr::select(starts_with('scDblFinder'))
  rownames(meta_scdblfinder) <- sce@colData@rownames
  cat("scDblFinder columns:", colnames(meta_scdblfinder), "\n")
  # reproducibility check: the same seed must give the same calls
  sce2 <- as.SingleCellExperiment(seu); set.seed(SEED_DBL); sce2 <- scDblFinder(sce2)
  same_class <- identical(as.character(sce2$scDblFinder.class), as.character(sce$scDblFinder.class))
  max_score_diff <- max(abs(sce2$scDblFinder.score - sce$scDblFinder.score))
  cat("re-run with the same seed: identical class:", same_class, "; max |score diff|:", max_score_diff, "\n")
  rm(sce2)
  m <- match(cells$cell, rownames(meta_scdblfinder))
  cells$nCount_RNA_noMT   <- NA; cells$nFeature_RNA_noMT <- NA
  cells$nCount_RNA_noMT[!is.na(m)]   <- Matrix::colSums(cf)[cells$cell[!is.na(m)]]
  cells$nFeature_RNA_noMT[!is.na(m)] <- Matrix::colSums(cf > 0)[cells$cell[!is.na(m)]]
  cells$scDblFinder_score <- meta_scdblfinder$scDblFinder.score[m]
  cells$scDblFinder_class <- as.character(meta_scdblfinder$scDblFinder.class)[m]
  sing <- rownames(meta_scdblfinder)[meta_scdblfinder$scDblFinder.class == "singlet"]
  out$cells <- cells; out$cts.filter <- cf; out$cts.sub <- cf[, sing, drop = FALSE]; out$scDblFinder <- meta_scdblfinder
  out$dbl_rerun_identical <- same_class; out$dbl_rerun_max_score_diff <- max_score_diff
  cat(sprintf("QC-pass %d | doublets %d (%.2f%%) | singlets %d\n", ncol(cf), sum(cells$scDblFinder_class == "doublet", na.rm = TRUE),
              100 * mean(meta_scdblfinder$scDblFinder.class == "doublet"), length(sing)))
} else {
  cat("LIBRARY FAILS the library rule: no doublet detection, no clustering.\n")
  cells$nCount_RNA_noMT <- NA; cells$nFeature_RNA_noMT <- NA; cells$scDblFinder_score <- NA; cells$scDblFinder_class <- NA
  out$cells <- cells
}
safe_save_rds(out, file.path(SCR, LIB, paste0(LIB, "_step01_qc_doublets.rds")))
sessionInfo()
