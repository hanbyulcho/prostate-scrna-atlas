# Code release header: estimates the ambient-RNA contamination fraction (rho) of one GEM well with SoupX 1.6.2. The counts are not corrected.
# Input: one Cell Ranger output folder. Settings: load10X(), then autoEstCont(doPlot = FALSE) at its defaults, after set.seed(<seed argument>).
# A well where SoupX fails is recorded as "not estimable", and nothing is retried. Output: one result row and the fit table per well.
# A script not included assembled the rows into qc_summary_sample_level.csv. Methods: Quality control.
# SoupX 1.6.2 contamination-fraction estimate (rho) for one Cell Ranger outs folder (one GEM well).
# Disclosure only; rho only; no adjustCounts; "not estimable" where SoupX fails.
# Method = SoupX vignette defaults: sc <- load10X(<outs>); sc <- autoEstCont(sc, doPlot = FALSE).
# No parameter is changed and nothing is retried.
#
# Usage: Rscript --no-environ 01_soupx_unit.R <unit_id> <outs_dir> <mode> <multi_sample> <seed> <result_dir>
#   mode = outs          : load10X(<outs_dir>)  (the main method)
#   mode = staged_multi  : Cell Ranger multi layout only. A folder of three symlinks is made in TMPDIR
#                          (raw_feature_bc_matrix -> multi/count/raw_feature_bc_matrix,
#                           filtered_feature_bc_matrix -> per_sample_outs/<s>/count/sample_filtered_feature_bc_matrix,
#                           analysis -> per_sample_outs/<s>/count/analysis) and load10X() is called on it.
#                          Secondary analysis, reported separately (not the main method).
#   multi_sample = per_sample_outs sample name (mode staged_multi), "-" otherwise.
# Writes <result_dir>/<unit_id>__<mode>__seed<seed>.csv (one row) and ..._fit_dd.csv.gz (sc$fit$dd + cluster).

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 6)
unit_id <- args[1]; outs_dir <- args[2]; mode <- args[3]; multi_sample <- args[4]
seed <- as.integer(args[5]); result_dir <- args[6]
stopifnot(mode %in% c("outs", "staged_multi"), !is.na(seed), dir.exists(outs_dir), dir.exists(result_dir))

R <- "<PROJECT_DIR>"
.libPaths(c(file.path(R, "redo_env_v1/Rlib_R433_overlay"), "<R_LIBRARY>", .Library))
suppressPackageStartupMessages(library(SoupX))
if (packageVersion("SoupX") != "1.6.2") stop("SoupX version is ", as.character(packageVersion("SoupX")), ", not 1.6.2")

script_path <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
script_md5 <- unname(tools::md5sum(script_path))
t0 <- Sys.time()
cat(sprintf("[%s] unit %s mode %s seed %d\n  outs: %s\n  script: %s (md5 %s)\n  TMPDIR: %s\n",
            format(t0), unit_id, mode, seed, outs_dir, script_path, script_md5, tempdir()))

# ---- load directory -------------------------------------------------------
load_dir <- outs_dir
if (mode == "staged_multi") {
  s_dir <- file.path(outs_dir, "per_sample_outs", multi_sample, "count")
  links <- c(raw_feature_bc_matrix = file.path(outs_dir, "multi", "count", "raw_feature_bc_matrix"),
             filtered_feature_bc_matrix = file.path(s_dir, "sample_filtered_feature_bc_matrix"),
             analysis = file.path(s_dir, "analysis"))
  stopifnot(all(dir.exists(links)))
  tmp_root <- Sys.getenv("TMPDIR")
  stopifnot(startsWith(tmp_root, "<SCRATCH_DIR>/soupx_tmp"))
  load_dir <- file.path(tmp_root, "staged_multi", unit_id)
  if (dir.exists(load_dir)) unlink(load_dir, recursive = TRUE)
  dir.create(load_dir, recursive = TRUE)
  for (nm in names(links)) stopifnot(file.symlink(links[[nm]], file.path(load_dir, nm)))
  cat("  staged load dir:", load_dir, "\n"); print(links)
}

# ---- run with condition capture --------------------------------------------
msgs <- character(0); warns <- character(0)
run_stage <- function(stage, expr) {
  withCallingHandlers(expr,
    message = function(m) { msgs <<- c(msgs, paste0(stage, ": ", trimws(conditionMessage(m))))
                            cat("  [message]", stage, trimws(conditionMessage(m)), "\n"); invokeRestart("muffleMessage") },
    warning = function(w) { warns <<- c(warns, paste0(stage, ": ", trimws(conditionMessage(w))))
                            cat("  [warning]", stage, trimws(conditionMessage(w)), "\n"); invokeRestart("muffleWarning") })
}

set.seed(seed)
rng_kind <- paste(RNGkind(), collapse = "/")
rs0 <- .Random.seed
status <- "estimated"; err_stage <- ""; err_text <- ""
sc <- tryCatch(run_stage("load10X", load10X(load_dir)),
               error = function(e) { status <<- "not estimable"; err_stage <<- "load10X"
                                     err_text <<- trimws(conditionMessage(e)); NULL })
rs1 <- .Random.seed
sc_loaded <- sc
if (!is.null(sc)) {
  sc <- tryCatch(run_stage("autoEstCont", autoEstCont(sc, doPlot = FALSE)),
                 error = function(e) { status <<- "not estimable"; err_stage <<- "autoEstCont"
                                       err_text <<- trimws(conditionMessage(e)); NULL })
}
rs2 <- .Random.seed
if (nzchar(err_text)) cat("  [error]", err_stage, err_text, "\n")

# ---- record ----------------------------------------------------------------
na <- NA
row <- list(unit_id = unit_id, mode = mode, outs_dir = outs_dir, load_dir = load_dir, seed = seed,
            rng_kind = rng_kind, status = status, error_stage = err_stage, error_text = err_text,
            warnings = paste(warns, collapse = " || "), messages = paste(msgs, collapse = " || "),
            rho = na, rho_fit_rhoEst = na, rho_metaData_unique = na, rho_fit_equals_metaData = na,
            rho_fwhm_low = na, rho_fwhm_high = na, prior_rho = na, prior_rho_sd = na,
            n_cells = na, n_genes = na, n_droplets_raw = na, n_droplets_soup_0_100 = na,
            n_clusters = na, cluster_source = na,
            n_markers_tfidf_pass = na, n_markers_soup_filter_pass = na, n_marker_genes_used = na,
            n_estimates_used = na, n_genes_with_estimate = na, n_clusters_with_estimate = na,
            rng_unchanged_by_load10X = identical(rs0, rs1), rng_unchanged_by_autoEstCont = identical(rs1, rs2),
            soupx_version = as.character(packageVersion("SoupX")), seurat_version = as.character(packageVersion("Seurat")),
            r_version = R.version.string, script_md5 = script_md5,
            started = format(t0, "%Y-%m-%d %H:%M:%S %Z"), elapsed_s = na)

if (!is.null(sc_loaded)) {
  row$n_cells <- ncol(sc_loaded$toc); row$n_genes <- nrow(sc_loaded$toc)
  row$n_droplets_raw <- length(sc_loaded$nDropUMIs)
  row$n_droplets_soup_0_100 <- sum(sc_loaded$nDropUMIs > 0 & sc_loaded$nDropUMIs < 100)  # estimateSoup default soupRange
  if ("clusters" %in% colnames(sc_loaded$metaData)) {
    row$n_clusters <- length(unique(sc_loaded$metaData$clusters))
    row$cluster_source <- "Cell Ranger graphclust (analysis/clustering/gene_expression_graphclust/clusters.csv), loaded by load10X"
  } else row$cluster_source <- "none loaded"
}
dd_out <- NULL
if (!is.null(sc)) {
  f <- sc$fit
  # autoEstCont can return without an error but with no estimate: with 0 usable gene/cluster estimates the
  # posterior is NaN everywhere, which.max() is empty, rhoEst = numeric(0) and no rho is set (SoupX 1.6.2 source).
  has_est <- length(f$rhoEst) == 1 && is.finite(f$rhoEst)
  rmeta <- unique(sc$metaData$rho)
  row$rho_metaData_unique <- if (length(rmeta)) paste(rmeta, collapse = ";") else NA
  row$prior_rho <- f$priorRho; row$prior_rho_sd <- f$priorRhoStdDev
  if (has_est) {
    row$rho_fit_rhoEst <- f$rhoEst
    row$rho_fit_equals_metaData <- length(rmeta) == 1 && identical(as.numeric(rmeta), as.numeric(f$rhoEst))
    if (!isTRUE(row$rho_fit_equals_metaData)) stop("sc$fit$rhoEst and sc$metaData$rho differ")
    row$rho <- f$rhoEst
    row$rho_fwhm_low <- f$rhoFWHM[1]; row$rho_fwhm_high <- f$rhoFWHM[2]
  } else {
    status <- "not estimable"; err_stage <- "autoEstCont (no estimate returned)"
    err_text <- sprintf("autoEstCont returned without an error, but sc$fit$rhoEst is empty (length %d) and no rho was set; %d independent estimates of rho were usable",
                        length(f$rhoEst), sum(f$dd$useEst))
    row$status <- status; row$error_stage <- err_stage; row$error_text <- err_text
    cat("  [no estimate]", err_text, "\n")
  }
  row$n_markers_tfidf_pass <- nrow(f$markersUsed)
  m <- regmatches(row$messages, regexec("(\\d+) genes passed tf-idf cut-off and (\\d+) soup quantile filter\\.\\s+Taking the top (\\d+)", row$messages))[[1]]
  if (length(m) == 4) {
    if (as.integer(m[2]) != row$n_markers_tfidf_pass) stop("tf-idf count in message differs from nrow(sc$fit$markersUsed)")
    row$n_markers_soup_filter_pass <- as.integer(m[3])
  }
  dd <- f$dd
  genes_used <- unique(dd$gene)
  cl_names <- names(split(rownames(sc$metaData), sc$metaData$clusters))  # column order used inside autoEstCont
  if (nrow(dd) != length(genes_used) * length(cl_names) || !identical(dd$gene, rep(genes_used, length(cl_names))))
    stop("unexpected layout of sc$fit$dd")
  dd$cluster <- rep(cl_names, each = length(genes_used))
  if (length(m) == 4 && as.integer(m[4]) != length(genes_used)) stop("'Taking the top' count differs from genes in sc$fit$dd")
  row$n_marker_genes_used <- length(genes_used)
  row$n_estimates_used <- sum(dd$useEst)
  row$n_genes_with_estimate <- length(unique(dd$gene[dd$useEst]))
  row$n_clusters_with_estimate <- length(unique(dd$cluster[dd$useEst]))
  m2 <- regmatches(row$messages, regexec("Using (\\d+) independent estimates of rho", row$messages))[[1]]
  if (length(m2) == 2 && as.integer(m2[2]) != row$n_estimates_used) stop("estimate count in message differs from sum(dd$useEst)")
  dd_out <- dd
}
row$elapsed_s <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)

# ---- safe writes: temp file -> read back -> check rows -> rename ------------
stem <- file.path(result_dir, sprintf("%s__%s__seed%d", unit_id, mode, seed))
safe_csv <- function(df, path, gz = FALSE) {
  tmp <- paste0(path, ".tmp")
  con <- if (gz) gzfile(tmp, "w") else file(tmp, "w")
  write.csv(df, con, row.names = FALSE); close(con)
  back <- read.csv(tmp, check.names = FALSE)
  if (nrow(back) != nrow(df) || ncol(back) != ncol(df)) stop("read-back mismatch for ", path)
  if (!file.rename(tmp, path)) stop("rename failed for ", path)
}
if (!is.null(dd_out)) safe_csv(dd_out, paste0(stem, "_fit_dd.csv.gz"), gz = TRUE)
safe_csv(as.data.frame(row, stringsAsFactors = FALSE), paste0(stem, ".csv"))
cat(sprintf("  status %s; rho %s; fwhm %s-%s; clusters %s; genes used %s; estimates %s; %.1f s\n",
            status, row$rho, row$rho_fwhm_low, row$rho_fwhm_high, row$n_clusters, row$n_marker_genes_used,
            row$n_estimates_used, row$elapsed_s))
if (mode == "staged_multi") unlink(load_dir, recursive = TRUE)
cat("---- sessionInfo ----\n"); print(sessionInfo())
