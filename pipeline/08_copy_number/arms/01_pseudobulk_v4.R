# Code release header: sums raw counts into pseudobulks of luminal units and of reference cell groups. Methods: Copy-number inference.
# Input: a cell table of kept cells and labels, its count-source manifest, each sample's count matrix and the cell metadata of the first atlas version.
# Output: pb_counts_v4.rds, pb_meta_v4.csv, halves_index_v4.csv and halves/, which 02_cnv_v4.R reads.
# A unit is the kept luminal cells of one cluster in a sample, if at least 50. Primary (D4) and comparison (D3) units use two clusterings.
# Only D4 units get 20 random half-splits. Reference groups: each sample's luminal cells and its non-epithelial cell types of 50 or more cells.
# 01_pseudobulk_v4.R
# Raw-count pseudobulks, all eligible units. Labels and kept cells: the stage-1 cell table (md5 200b2258783706a2042ff6f3f325a9e5).
# Counts: each cell's source h5 by library_id + gem_well in count_source_manifest_v12_all_v6.csv, barcode = source_barcode.
# D4 units (primary): library x cluster_source x cluster_res3, kept Luminal cells, >= 50 cells.
# D3 units (comparison): library x cluster_own_res3 of the label-check tables, kept Luminal cells, >= 50 cells.
# Half-splits: 20 per D4 unit. Tumor-type share: metacomb etsgroup in ERG, ETV, nonETS, AR, NEPC, DN; cells without a row there count as not tumor-type.
.libPaths(c("<PROJECT_DIR>/redo_env_v1/Rlib_R433_overlay",
            "<R_LIBRARY>", .Library))
suppressPackageStartupMessages({ library(data.table); library(Matrix); library(rhdf5) })
t0 <- Sys.time(); args <- commandArgs(TRUE); TESTMODE <- length(args) > 0 && args[1] == "test"
R   <- "<PROJECT_DIR>"
B   <- "<BASE_DIR>"
W   <- file.path("<work folder: the script writes the pseudobulks and half-splits under it>", if (TESTMODE) "work_test" else "work"); dir.create(file.path(W, "halves"), showWarnings = FALSE, recursive = TRUE)
CLT <- "<input folder: holds <library>_cell_clusters.csv.gz, the per-library cluster tables of the label check>"
MIN_UNIT <- 50L; NSPLIT <- 20L
EPI <- c("Luminal","Basal","Club","Hillock","Ionocyte"); TUM <- c("ERG","ETV","nonETS","AR","NEPC","DN")
SRC <- function(x) fifelse(grepl("^first", x), "fs", fifelse(grepl("^recall", x), "rc", fifelse(grepl("^new_lib", x), "nl", fifelse(grepl("^ts_per", x), "ts", NA_character_))))

C <- fread(file.path(R, "release_build_v1/outputs/dryrun_v12/cell_table_stage1_dryrun_v12_v6.csv.gz"),
           select = c("cell_id","library_id","gem_well","source_barcode","cell_id_ts_per_channel_v1","count_source_md5","first_submission_cell_id","cluster_res3","cluster_source","release_status","cell_type"))[release_status == "kept"]
C[, source := SRC(cluster_source)]; stopifnot(!anyNA(C$source), !anyNA(C$cluster_res3))
Mf <- fread(file.path(R, "release_build_v1/outputs/dryrun_v12/count_source_manifest_v12_all_v6.csv"), select = c("library_id","gem_well","path","md5"))
stopifnot(!anyDuplicated(Mf[, .(library_id, gem_well)]))
chk <- merge(unique(C[, .(library_id, gem_well, count_source_md5)]), Mf, by = c("library_id","gem_well"), all.x = TRUE); stopifnot(!anyNA(chk$path), all(chk$count_source_md5 == chk$md5))
hdr <- names(fread(file.path(B, "17datasets/metacomb_17data_final.csv"), nrows = 0)); ie <- which(hdr == "etsgroup"); stopifnot(length(ie) == 1)
MC <- fread(file.path(B, "17datasets/metacomb_17data_final.csv"), select = c(1L, ie)); setnames(MC, c("fs_id","etsgroup"))
C[, etsgroup := MC$etsgroup[match(first_submission_cell_id, MC$fs_id)]]; C[, tumour_type := etsgroup %in% TUM]
cat("luminal kept cells:", C[cell_type == "Luminal", .N], "with a metacomb row:", C[cell_type == "Luminal" & !is.na(etsgroup), .N], "tumor-type:", C[cell_type == "Luminal", sum(tumour_type)], "\n")
LIBS <- sort(unique(C[cell_type == "Luminal"]$library_id)); if (TESTMODE) LIBS <- c("L009","L138")

PB <- list(); META <- list(); CHECK <- list(); HI <- list(); LC <- list(); gene_ids <- NULL
read10x <- function(p) {
  ft <- as.vector(h5read(p, "/matrix/features/feature_type")); id <- as.vector(h5read(p, "/matrix/features/id"))
  bc <- as.vector(h5read(p, "/matrix/barcodes")); sh <- as.integer(h5read(p, "/matrix/shape"))
  ip <- as.integer(h5read(p, "/matrix/indptr")); ix <- as.integer(h5read(p, "/matrix/indices")); x <- as.numeric(h5read(p, "/matrix/data"))
  stopifnot(!anyNA(ip), length(ip) == sh[2] + 1L)
  m <- Matrix::sparseMatrix(i = ix + 1L, p = ip, x = x, dims = sh); rm(ix, x); invisible(gc())
  dimnames(m) <- list(id, bc); m[which(ft == "Gene Expression"), , drop = FALSE] }
for (L in LIBS) {
  K <- C[library_id == L]; wells <- unique(K$gem_well)
  ms <- lapply(wells, function(w) { x <- read10x(Mf[library_id == L & gem_well == w]$path); colnames(x) <- paste(w, colnames(x), sep = "|"); x })
  if (is.null(gene_ids)) gene_ids <- rownames(ms[[1]]); for (x in ms) stopifnot(identical(rownames(x), gene_ids))
  m <- if (length(ms) == 1) ms[[1]] else do.call(cbind, ms); rm(ms)
  K[, key := paste(gem_well, source_barcode, sep = "|")]; stopifnot(!anyDuplicated(K$key))
  hit <- K$key %in% colnames(m)
  lum <- K[cell_type == "Luminal"]
  CHECK[[L]] <- data.table(library = L, wells = paste(wells, collapse = ";"), n_kept = nrow(K), n_found = sum(hit), n_lum = nrow(lum)); stopifnot(all(hit))
  ## D3: label-check per-library clustering (cell id = per-channel id for TS)
  ft3 <- file.path(CLT, paste0(L, "_cell_clusters.csv.gz"))
  if (file.exists(ft3)) { cl <- fread(ft3, select = c("cell","cluster_own_res3"))
    lum[, cid3 := fifelse(!is.na(cell_id_ts_per_channel_v1) & cell_id_ts_per_channel_v1 != "", cell_id_ts_per_channel_v1, cell_id)]
    lum[, c3 := as.character(cl$cluster_own_res3[match(cid3, cl$cell)])] } else lum[, c3 := NA_character_]
  lum[, c4 := paste0(source, "|c", cluster_res3)]
  grp <- list()
  add <- function(id, kind, source, cluster, celltype, d) grp[[id]] <<- list(kind = kind, source = source, cluster = cluster, celltype = celltype, key = d$key, ntum = sum(d$tumour_type), nolab = sum(is.na(d$etsgroup)))
  add(paste0(L, "|lib_lum"), "lib_luminal", NA_character_, NA_character_, "Luminal", lum)
  for (k in sort(unique(lum$c4))) { d <- lum[c4 == k]; if (nrow(d) >= MIN_UNIT) add(paste0(L, "|", k), "unit_d4", d$source[1], as.character(d$cluster_res3[1]), "Luminal", d) }
  for (k in sort(unique(na.omit(lum$c3)))) { d <- lum[c3 == k]; if (nrow(d) >= MIN_UNIT) add(paste0(L, "|o", k), "unit_d3", "own", k, "Luminal", d) }
  ne <- K[!cell_type %in% EPI]
  for (ct in sort(unique(ne$cell_type))) { d <- ne[cell_type == ct]; if (nrow(d) >= MIN_UNIT) add(paste0(L, "|ne_", gsub(" ", "_", ct)), "nonepi_celltype", NA_character_, NA_character_, ct, d) }
  for (g in names(grp)) { b <- grp[[g]]$key; PB[[g]] <- as.numeric(Matrix::rowSums(m[, b, drop = FALSE]))
    META[[g]] <- data.table(pb_id = g, library = L, kind = grp[[g]]$kind, source = grp[[g]]$source, cluster = grp[[g]]$cluster, celltype = grp[[g]]$celltype,
                            n_cells = length(b), total_umi = sum(PB[[g]]), n_tumour_type = grp[[g]]$ntum, n_no_fs_row = grp[[g]]$nolab) }
  lum[, d4_unit := fifelse(paste0(L, "|", c4) %in% names(grp), paste0(L, "|", c4), NA_character_)]; lum[, d3_unit := fifelse(!is.na(c3) & paste0(L, "|o", c3) %in% names(grp), paste0(L, "|o", c3), NA_character_)]
  LC[[L]] <- lum[, .(cell_id, library = L, d4_unit, d3_unit, c3_missing = is.na(c3))]
  ## half-splits for D4 units (seed = split x 1e5 + unit index within the library; project choice)
  units <- names(grp)[vapply(grp, function(z) z$kind == "unit_d4", logical(1))]
  if (length(units)) { H <- list()
    for (ui in seq_along(units)) { b <- grp[[units[ui]]]$key; n <- length(b); ii <- integer(0); jj <- integer(0)
      for (s in seq_len(NSPLIT)) { set.seed(s * 1e5 + ui); h1 <- sort(sample.int(n, floor(n / 2))); h2 <- setdiff(seq_len(n), h1)
        ii <- c(ii, h1, h2); jj <- c(jj, rep(2L * s - 1L, length(h1)), rep(2L * s, length(h2))) }
      I <- Matrix::sparseMatrix(i = ii, j = jj, x = 1, dims = c(n, 2L * NSPLIT))
      hm <- as.matrix(m[, b, drop = FALSE] %*% I); storage.mode(hm) <- "integer"
      colnames(hm) <- sprintf("%s|s%02d|h%d", units[ui], rep(seq_len(NSPLIT), each = 2), rep(1:2, NSPLIT)); H[[ui]] <- hm
      HI[[units[ui]]] <- data.table(pb_id = units[ui], n_cells = n, n_half1 = floor(n / 2)) }
    H <- do.call(cbind, H); rownames(H) <- gene_ids
    f <- file.path(W, "halves", paste0(L, ".rds")); saveRDS(H, paste0(f, ".tmp")); system2("sync", paste0(f, ".tmp")); stopifnot(file.size(paste0(f, ".tmp")) > 0); file.rename(paste0(f, ".tmp"), f) }
  cat(sprintf("%s %s wells %d kept %d lum %d d4 %d d3 %d ne %d c3_missing %d\n", format(Sys.time(), "%T"), L, length(wells), nrow(K), nrow(lum), length(units),
              sum(vapply(grp, function(z) z$kind == "unit_d3", logical(1))), sum(vapply(grp, function(z) z$kind == "nonepi_celltype", logical(1))), sum(is.na(lum$c3))))
  rm(m); invisible(gc())
}
X <- do.call(cbind, PB); rownames(X) <- gene_ids; META <- rbindlist(META); CHECK <- rbindlist(CHECK); HI <- rbindlist(HI); LC <- rbindlist(LC)
stopifnot(identical(META$pb_id, colnames(X)))
sw <- function(obj, f, writer, nr = NULL) { tmp <- paste0(f, ".tmp"); writer(obj, tmp); system2("sync", tmp); stopifnot(file.size(tmp) > 0); if (!is.null(nr)) stopifnot(nrow(fread(tmp)) == nr); file.rename(tmp, f) }
sw(X, file.path(W, "pb_counts_v4.rds"), function(o, f) saveRDS(o, f))
fw <- function(o, f) fwrite(o, f); fg <- function(o, f) fwrite(o, f, compress = "gzip")
sw(META, file.path(W, "pb_meta_v4.csv"), fw, nrow(META)); sw(CHECK, file.path(W, "pb_cell_match_v4.csv"), fw, nrow(CHECK)); sw(HI, file.path(W, "halves_index_v4.csv"), fw, nrow(HI))
sw(LC, file.path(W, "luminal_cells_units_v4.csv.gz"), fg, nrow(LC))
cat("pseudobulks:", ncol(X), "\n"); print(META[, .N, by = kind]); cat("luminal cells without a D3 cluster:", sum(LC$c3_missing), "\n")
cat("elapsed_s", as.numeric(difftime(Sys.time(), t0, units = "secs")), "\n")
