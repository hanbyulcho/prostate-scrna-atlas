# Code release header: removes barcode-swapped molecules from all samples of one study with DropletUtils swappedDrops (min.frac = 0.8).
# Input: each sample's molecule_info.h5 and filtered matrix, listed in a sample table (placeholder folder). Methods: Raw sequencing processing.
# Cell QC on the cleaned counts: at least 250 genes and 500 UMIs, and under 20% mitochondrial UMIs. A sample needs at least 100 passing cells.
# Output: cleaned counts of the Cell Ranger-called barcodes and a removal reason per cell. The released counts of these studies come from 03_recall_qc_classify.R.
# swappedDrops (DropletUtils; Griffiths et al. 2018, Nat Commun 9:2667, doi:10.1038/s41467-018-05083-x) on ALL
# molecule_info.h5 of one study, then the fixed cell QC
# (nFeature >= 250 & nCount >= 500 & percent.mt < 20; metrics on all genes; MT = the
# 13 'MT-' genes) on the cleaned counts, then the library rule (>= 100 Cell Ranger-called barcodes passing that
# QC, before doublets). emptyDrops later re-called cells after the same swap removal on all
# barcodes; the released counts are those written by pipeline/02_qc/recall/03_recall_qc_classify.R. A molecule is kept in the library holding at least 80% of its reads.
#
# Checks and fixes (labels used below): F-A the written h5 files carry the
#   Cell Ranger genome, the library's chemistry and library id (DropletUtils 1.22.0 write10xCounts() does not pass
#   chemistry/library.ids to the file, so the root attributes are rewritten with rhdf5 and asserted); F-B entry-wise check
#   of the counts before swap removal against the Cell Ranger matrix, and the re-read h5 must equal the in-memory matrix entry by
#   entry; the optional comparison table (CMP argument) makes the script exit with status 2 when any mismatch is found (outputs are written first).
# Further checks and fixes:
#   F1 cleaned counts of the called barcodes written per library (10x HDF5 v3), re-read and checked, then renamed;
#   F2 counts checked against the Cell Ranger filtered matrix; every atlas or raw-object cell of the first atlas version must pass the QC on the Cell Ranger counts; the
#      QC pass count on the Cell Ranger counts must equal qc_Cells_after of the first atlas version where recorded; whole-study grouping enforced;
#   F4 removal reason per cell ('qc_after_swap:<criteria>' / 'library_rule'), no double spaces; library_rule_pass per cell;
#   F5 data.table set() (no selfref warning); script path + md5 and molecule_info paths logged.
#
# Usage: Rscript --no-environ 04_swappeddrops_study.R STUDY LIB1,LIB2,... OUTDIR [CMP_CELLS_CSV]
REV <- "<PROJECT_DIR>"
.libPaths(c(file.path(REV, "redo_env_v1/Rlib_R433_overlay"), "<R_LIBRARY>", .Library))
suppressPackageStartupMessages({library(DropletUtils); library(Matrix); library(data.table); library(rhdf5); library(jsonlite)})

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 3)
STUDY <- args[1]; LIBS <- strsplit(args[2], ",")[[1]]; OUT <- args[3]; CMP <- if (length(args) >= 4) args[4] else NA
E2 <- "<input folder: holds tables/library_table.csv (one row per library), tmp/obs_integrated.csv.gz and tmp/obs_raw.csv.gz (cell IDs of the first-version integrated and raw objects)>"
SELF <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
cat("command:", paste(commandArgs(), collapse = " "), "\n")
cat("script:", SELF, "md5", tools::md5sum(SELF), "\n")
cat("libPaths:", paste(.libPaths(), collapse = " | "), "\n")
cat("TMPDIR:", Sys.getenv("TMPDIR"), "\n")
cat("start:", format(Sys.time()), "\n")
dir.create(file.path(OUT, "counts"), showWarnings = FALSE, recursive = TRUE)

safe_fwrite <- function(x, path) {
  tmp <- paste0(path, ".tmp", if (grepl("\\.gz$", path)) ".gz" else "")
  fwrite(x, tmp)
  stopifnot(file.exists(tmp), file.size(tmp) > 0, nrow(fread(tmp)) == nrow(x))
  stopifnot(file.rename(tmp, path))
}

## inputs and grouping (whole study, including excluded/pending libraries)
lt_all <- fread(file.path(E2, "tables/library_table.csv"))
stopifnot(!anyDuplicated(LIBS), all(LIBS %in% lt_all$ID))
PAPER <- unique(lt_all[ID %in% LIBS, Paper_ID])
stopifnot(length(PAPER) == 1)
if (!setequal(LIBS, lt_all[Paper_ID == PAPER, ID])) stop("grouping must be the whole study ", PAPER, ": ", paste(lt_all[Paper_ID == PAPER, ID], collapse = ","))
lt <- lt_all[match(LIBS, ID)]
files <- setNames(lt$molecule_info, LIBS); filt <- setNames(lt$filtered_h5, LIBS)
stopifnot(all(file.exists(files)), all(file.exists(filt)))
for (L in LIBS) cat("input", L, files[[L]], "|", filt[[L]], "\n")

## per-file asserts
feat_id <- NULL; chem <- setNames(character(length(LIBS)), LIBS)
for (L in LIBS) {
  li <- fromJSON(h5read(files[[L]], "library_info")); mj <- fromJSON(h5read(files[[L]], "metrics_json"))
  ft <- as.character(h5read(files[[L]], "features/feature_type")); id <- as.character(h5read(files[[L]], "features/id"))
  if (is.null(feat_id)) { feat_id <- id; feat_name <- as.character(h5read(files[[L]], "features/name")) }
  ok <- c(one_library = nrow(li) == 1, gene_expression = all(li$library_type == "Gene Expression") && all(ft == "Gene Expression"),
          one_gem_group = length(unique(li$gem_group)) == 1, same_features = identical(id, feat_id),
          include_introns = isTRUE(mj$analysis_parameters$include_introns))
  chem[[L]] <- mj$chemistry_description
  cat(L, lt[ID == L, run], paste(names(ok), ok, sep = "=", collapse = " "), chem[[L]], "\n")
  if (!all(ok)) stop("input assert failed for ", L)
}
h5closeAll()
mt_id <- feat_id[grepl("^MT-", feat_name)]
cat("MT genes:", length(mt_id), paste(feat_name[grepl("^MT-", feat_name)], collapse = ","), "\n")
stopifnot(length(mt_id) == 13)

atlas <- fread(file.path(E2, "tmp/obs_integrated.csv.gz"), select = 1)[[1]]
rawset <- fread(file.path(E2, "tmp/obs_raw.csv.gz"), select = 1)[[1]]

## swappedDrops
t0 <- Sys.time()
s <- swappedDrops(files, min.frac = 0.8, get.swapped = TRUE, get.diagnostics = FALSE)
cat("swappedDrops (DropletUtils", as.character(packageVersion("DropletUtils")), ") done in", format(Sys.time() - t0), "\n")
stopifnot(identical(names(s$cleaned), LIBS), identical(names(s$swapped), LIBS))

cells <- list(); libsum <- list()
for (L in LIBS) {
  cl <- s$cleaned[[L]]; sw <- s$swapped[[L]]
  stopifnot(identical(dim(cl), dim(sw)), identical(colnames(cl), colnames(sw)), all(rownames(cl) == feat_id), nrow(cl) == length(feat_id))
  ## Cell Ranger filtered matrix: called barcodes and their counts (F2)
  called <- as.character(h5read(filt[[L]], "matrix/barcodes"))
  cr_ids <- as.character(h5read(filt[[L]], "matrix/features/id"))
  cr_data <- as.numeric(h5read(filt[[L]], "matrix/data")); cr_ptr <- as.numeric(h5read(filt[[L]], "matrix/indptr"))
  cr_idx <- as.integer(h5read(filt[[L]], "matrix/indices"))
  gen <- unique(as.character(h5read(filt[[L]], "matrix/features/genome"))); h5closeAll(); stopifnot(length(gen) == 1)
  stopifnot(identical(cr_ids, feat_id), length(cr_ptr) == length(called) + 1)
  colsum_csc <- function(x, ptr) { cs <- c(0, cumsum(x)); cs[ptr[-1] + 1] - cs[ptr[-length(ptr)] + 1] }
  cr_umi <- colsum_csc(cr_data, cr_ptr); cr_genes <- as.integer(diff(cr_ptr))
  mtrow <- (cr_idx + 1L) %in% which(feat_id %in% mt_id)
  cr_mt <- colsum_csc(ifelse(mtrow, cr_data, 0), cr_ptr)
  called_bc <- sub("-1$", "", called)
  if (L == LIBS[1]) cat("cleaned colnames example:", head(colnames(cl), 2), "| called example:", head(called, 2), "\n")
  stopifnot(length(called) > 0, all(grepl("-1$", called)), !anyDuplicated(called_bc), all(called_bc %in% colnames(cl)))
  cl <- cl[, called_bc, drop = FALSE]; sw <- sw[, called_bc, drop = FALSE]; og <- cl + sw
  mt <- rownames(cl) %in% mt_id
  d <- data.table(lib = L, barcode = called, cell = paste0(L, "_", called),
                  umi_orig = colSums(og), genes_orig = colSums(og > 0), mt_orig = colSums(og[mt, , drop = FALSE]),
                  umi_clean = colSums(cl), genes_clean = colSums(cl > 0), mt_clean = colSums(cl[mt, , drop = FALSE]),
                  umi_swapped = colSums(sw))
  ## F2: molecule-level counts before swap removal equal the Cell Ranger filtered matrix
  stopifnot(all(d$umi_orig == cr_umi), all(d$genes_orig == cr_genes), all(d$mt_orig == cr_mt))
  d[, pct_mt_orig := 100 * mt_orig / pmax(umi_orig, 1)][, pct_mt_clean := 100 * mt_clean / pmax(umi_clean, 1)]
  d[, pass_qc_orig := genes_orig >= 250 & umi_orig >= 500 & pct_mt_orig < 20]
  d[, pass_qc_clean := genes_clean >= 250 & umi_clean >= 500 & pct_mt_clean < 20]
  crit <- cbind(ifelse(d$genes_clean < 250, "genes<250", ""), ifelse(d$umi_clean < 500, "umi<500", ""), ifelse(d$pct_mt_clean >= 20, "mt>=20%", ""))
  d[, fail_reason_clean := apply(crit, 1, function(z) paste(z[z != ""], collapse = ";"))]
  d[, in_atlas := cell %in% atlas][, in_raw_object := cell %in% rawset]
  stopifnot(all(atlas[startsWith(atlas, paste0(L, "_"))] %in% d$cell), all(rawset[startsWith(rawset, paste0(L, "_"))] %in% d$cell))
  crm <- sparseMatrix(i = cr_idx + 1L, p = as.integer(cr_ptr), x = cr_data, dims = c(length(feat_id), length(called)))
  stopifnot(sum(abs(og - crm)) == 0)                      # F-B: cleaned + swapped equals the Cell Ranger matrix entry by entry
  rm(crm)
  ## F2: every cell of the first atlas version passes the QC on the Cell Ranger counts; the pass count equals its recorded number
  stopifnot(all(d$pass_qc_orig[d$in_atlas | d$in_raw_object]))
  qa <- lt[ID == L, qc_Cells_after]
  if (!is.na(qa)) stopifnot(sum(d$pass_qc_orig) == qa)
  n_pass <- sum(d$pass_qc_clean); rule <- n_pass >= 100
  d[, library_rule_pass := rule]
  d[, removal_reason := fifelse(!pass_qc_clean, paste0("qc_after_swap:", fail_reason_clean), fifelse(!library_rule_pass, "library_rule", ""))]
  ## F1: write the cleaned counts of the called barcodes (10x HDF5 v3), re-read, check, rename
  h5 <- file.path(OUT, "counts", paste0(L, "_swappeddrops_cleaned_filtered_feature_bc_matrix.h5")); tmp <- sub("\\.h5$", ".tmp.h5", h5)
  if (file.exists(tmp)) file.remove(tmp)
  colnames(cl) <- called
  write10xCounts(tmp, cl, barcodes = called, gene.id = feat_id, gene.symbol = feat_name, type = "HDF5", version = "3",
                 genome = gen, overwrite = TRUE)
  fid <- H5Fopen(tmp); gid <- H5Gopen(fid, "/")         # F-A: write10xCounts() does not pass chemistry / library ids to the file
  for (nm in c("chemistry_description", "library_ids")) if (H5Aexists(gid, nm)) H5Adelete(gid, nm)
  h5writeAttribute(chem[[L]], h5obj = gid, name = "chemistry_description", variableLengthString = TRUE, asScalar = TRUE, encoding = "UTF-8")
  h5writeAttribute(L, h5obj = gid, name = "library_ids", variableLengthString = TRUE, asScalar = TRUE, encoding = "UTF-8")
  H5Gclose(gid); H5Fclose(fid)
  rb <- as.character(h5read(tmp, "matrix/barcodes")); rp <- as.numeric(h5read(tmp, "matrix/indptr")); rd <- as.numeric(h5read(tmp, "matrix/data"))
  rid <- as.character(h5read(tmp, "matrix/features/id")); ri <- as.integer(h5read(tmp, "matrix/indices"))
  rnm <- as.character(h5read(tmp, "matrix/features/name")); rgen <- unique(as.character(h5read(tmp, "matrix/features/genome")))
  ra <- h5readAttributes(tmp, "/"); h5closeAll()
  stopifnot(identical(ri, cl@i), all(rd == cl@x), identical(rnm, feat_name), identical(rgen, gen),
            identical(as.character(ra$chemistry_description), chem[[L]]), identical(as.character(ra$library_ids), L))
  stopifnot(identical(rb, called), identical(rid, feat_id), length(rp) == length(rb) + 1,
            all(colsum_csc(rd, rp) == d$umi_clean), all(diff(rp) == d$genes_clean), sum(rd) == sum(d$umi_clean))
  stopifnot(file.rename(tmp, h5))
  libsum[[L]] <- data.table(study = STUDY, paper_id = PAPER, lib = L, run = lt[ID == L, run], exclude_meta = trimws(lt[ID == L, Exclude_meta]),
                            called = nrow(d), pass_qc_orig = sum(d$pass_qc_orig), qc_cells_after_recorded = qa, pass_qc_clean = n_pass,
                            library_rule_ge100 = rule, umi_swapped_total = sum(d$umi_swapped),
                            atlas_cells = sum(d$in_atlas), atlas_fail_qc_clean = sum(d$in_atlas & !d$pass_qc_clean),
                            atlas_removed_total = sum(d$in_atlas & d$removal_reason != ""),
                            raw_only_cells = sum(d$in_raw_object & !d$in_atlas),
                            raw_only_removed = sum(d$in_raw_object & !d$in_atlas & d$removal_reason != ""),
                            counts_h5 = h5, counts_h5_md5 = unname(tools::md5sum(h5)))
  cells[[L]] <- d
  cat(L, "called", nrow(d), "pass_clean", n_pass, "rule", rule, "atlas", sum(d$in_atlas), "atlas_removed", sum(d$in_atlas & d$removal_reason != ""), "\n")
}
cells <- rbindlist(cells); libsum <- rbindlist(libsum)
safe_fwrite(cells, file.path(OUT, paste0(STUDY, "_cells.csv.gz")))
safe_fwrite(libsum, file.path(OUT, paste0(STUDY, "_library_summary.csv")))
print(libsum[, !c("counts_h5", "counts_h5_md5")])
cat("atlas cells removed (QC after swap or library rule):", sum(cells$in_atlas & cells$removal_reason != ""), "of", sum(cells$in_atlas), "\n")

## per-cell comparison with the CMP table of an independent re-implementation (min_frac 0.8 columns)
if (!is.na(CMP)) {
  e <- fread(CMP, select = c("cell", "orig_umi", "orig_genes", "orig_mt", "umi08", "genes08", "mt08", "pass_qc08"))
  m <- merge(cells, e, by = "cell", all = TRUE)
  cmp <- data.table(study = STUDY, cells_v4 = nrow(cells), cells_14oct = nrow(e), in_both = sum(!is.na(m$umi_clean) & !is.na(m$umi08)),
                    only_v3 = sum(is.na(m$umi08)), only_14oct = sum(is.na(m$umi_clean)))
  b <- m[!is.na(umi_clean) & !is.na(umi08)]
  for (p in list(c("umi_orig", "orig_umi"), c("genes_orig", "orig_genes"), c("mt_orig", "orig_mt"), c("umi_clean", "umi08"),
                 c("genes_clean", "genes08"), c("mt_clean", "mt08"), c("pass_qc_clean", "pass_qc08")))
    set(cmp, j = paste0("mismatch_", p[1]), value = sum(b[[p[1]]] != b[[p[2]]]))
  set(cmp, j = "atlas_mismatch_pass", value = sum(b$in_atlas & b$pass_qc_clean != b$pass_qc08))
  safe_fwrite(cmp, file.path(OUT, paste0(STUDY, "_vs_14oct.csv")))
  print(t(cmp))
  diffs <- unlist(cmp[, c("only_v3", "only_14oct", grep("^mismatch_|^atlas_mismatch", names(cmp), value = TRUE)), with = FALSE])
  cat(if (all(diffs == 0)) "COMPARISON_OK: identical to the 14 Oct re-implementation\n" else "COMPARISON_DIFF: see *_vs_14oct.csv\n")
}
cat("end:", format(Sys.time()), "\n")
print(sessionInfo())
if (!is.na(CMP) && exists("diffs") && any(diffs != 0)) quit(status = 2)
