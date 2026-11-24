# Code release header: applies the emptyDrops calls to the swap-cleaned counts of all barcodes of one study, then the cell QC.
# Input: the all-barcode cleaned counts of step 1, and the emptyDrops calls of a separate check and of step 2 (none included).
#   It also reads the cell table of 04_swappeddrops_study.R and the cell metadata of the first atlas version.
#   It also reads F_UNITS (cells of sets A and B), F_NOGRP (raw-object-only Luminal_No_Group cells) and the doublet classes of the first atlas version.
#   For L071 it reads the cell table of a per-sample step not included (new_libraries_v1/L071/cells_L071.csv.gz).
# Cell QC: at least 250 genes and 500 UMIs, and under 20% mitochondrial UMIs. A sample needs at least 100 re-called cells that pass.
# Output: the re-called count matrix of each sample (10x HDF5) and a class for every Cell Ranger-called or re-called barcode. Methods: Raw sequencing processing, Quality control.
# 03_recall_qc_classify.R STUDY -- step 3 of the cell re-call: re-called set, fixed QC, library rule, pre-stated stop, barcode classes,
# re-called count matrices. Reads the all-barcode RDS of step 1 (several GB in all) -> run through a wrapper script (not included).
#   Calls: the DropletUtils-default calls (FDR <= 0.01) of a separate emptyDrops check, except libraries with uncalled Limited = TRUE barcodes, which take the
#     niters = 30000 re-run (step 2, mode niters). Step 2 mode verify must exist for one library and has asserted that
#     the default re-run equals that check's calls.
#   Fixed QC (as in 04_swappeddrops_study.R): genes >= 250 & UMI >= 500 & pct MT < 20 on all
#     genes of the cleaned counts; MT = the 13 'MT-' genes (= element 2 of <BASE_DIR>/basic.info.rds, asserted).
#   Library rule (counted on the re-called set): re-called barcodes passing the fixed QC >= 100.
#   Pre-stated stop: computed for every library of the group BEFORE any library is processed; if L161 or any of
#     L159, L160, L169, L170, L171 passes, the script writes the library-rule table and the stop report and halts (exit 5).
#   Barcode classes for every barcode that is Cell Ranger-called (call list of the first atlas version =
#     filtered_h5 of library_table.csv, = the swap step's cell rows) or re-called:
#     A  kept atlas cell (in_atlas & removal_reason == "") re-called, library passes  -> candidate to keep (label carried)
#     B  kept atlas cell not re-called                                                -> removal row "not re-called by emptyDrops after swap removal"
#     B2 kept atlas cell re-called, library fails the library rule on the re-called set     -> removal row "library_rule (re-called set)"
#     C  atlas cell already removed by the swap step (removal_reason != "")                  -> stays in the swap step's removal table (no new row)
#     D  Cell Ranger-called, not an atlas cell                                     -> keeps its fate in the first atlas version
#     E  newly called (re-called, not Cell Ranger-called)                          -> candidate new cell; QC, rule, doublets in step 4
#     F  L071 only: re-called barcode -> whole per-sample processing; G  L071 Cell Ranger-called, not re-called
# Outputs: <output folder: the folder OUT that 00_env.R sets>/<STUDY>/counts/<LIB>_swappeddrops_cleaned_recalled_feature_bc_matrix.h5 (cleaned counts of the re-called
#   barcodes, 10x HDF5 v3, written as 04_swappeddrops_study.R F1/F-A); tables/<STUDY>_03_*.
source("<PROJECT_DIR>/recall_v1/scripts/00_env.R")
suppressPackageStartupMessages({library(DropletUtils); library(Matrix); library(data.table); library(rhdf5); library(jsonlite)})
STUDY <- commandArgs(TRUE)[1]; stopifnot(STUDY %in% names(GROUPS)); PAPER <- GROUPS[[STUDY]]
script_header(); pkg_versions(); check_disk("start")
TMPD <- file.path(SCRB, STUDY); TB <- file.path(OUT, "tables"); OD <- file.path(OUT, STUDY); dir.create(file.path(OD, "counts"), recursive = TRUE, showWarnings = FALSE)
RUN4 <- file.path(REV, "swap_removal_v1/run_v4", STUDY)
EDC <- if (STUDY == "P09_Cheng") file.path(REV, "swap_removal_v1/emptydrops_check_cheng_v1/tables") else file.path(REV, "swap_removal_v1/emptydrops_check_other_groups_v1/tables")
F_BC   <- if (STUDY == "P09_Cheng") file.path(EDC, "02_emptydrops_barcodes_all_runs.csv.gz") else file.path(EDC, paste0(STUDY, "_02_emptydrops_barcodes.csv.gz"))
F_META <- if (STUDY == "P09_Cheng") file.path(EDC, "02_emptydrops_run_metadata.csv") else file.path(EDC, paste0(STUDY, "_02_emptydrops_run_metadata.csv"))
F_CELLS <- file.path(RUN4, paste0(STUDY, "_cells.csv.gz")); F_LT <- "<input: library table, one row per library (columns ID, Paper_ID, filtered_h5, molecule_info)>"
F_META1 <- file.path(C17, "metacomb_17data_final.csv"); F_UNITS <- "<input: cells of the curation units (columns lib, cell, set, unit_id; set A = audit units, set B = removal units)>"
F_NOGRP <- "<input: IDs of the raw-object-only Luminal_No_Group cells, one per line>"; F_DBL1 <- "<input: first-version scDblFinder class per cell (columns lib, cell, scDblFinder.class)>"
F_BASIC <- file.path(C, "basic.info.rds")
for (f in c(F_BC, F_META, F_CELLS, F_LT, F_META1, F_UNITS, F_NOGRP, F_DBL1, F_BASIC)) cat("input", f, "md5", md5(f), "\n")

lt <- fread(F_LT)[Paper_ID == PAPER]; LIBS <- lt$ID
man <- fread(file.path(TMPD, "allbarcodes_manifest.csv")); stopifnot(identical(man$lib, LIBS))
rv <- fread(F_BC); rmeta <- fread(F_META)[run == "cleaned_i"]; stopifnot(identical(rmeta$lib, LIBS), all(rmeta$seed == SEED_ED))
NITERS_LIBS <- rmeta[limited_not_called > 0, lib]                  # libraries with uncalled Limited = TRUE barcodes
cat("libraries with Limited = TRUE barcodes not called at the defaults:", if (length(NITERS_LIBS)) paste(NITERS_LIBS, collapse = ",") else "none", "\n")
VER <- list.files(TB, pattern = paste0("^", STUDY, "_02_ed_.*_verify_meta\\.csv$"))
stopifnot(length(VER) >= 1); cat("default re-run verification(s):", paste(VER, collapse = ", "), "\n")
for (L in NITERS_LIBS) stopifnot(file.exists(file.path(TB, paste0(STUDY, "_02_ed_", L, "_niters.csv.gz"))))
cells4 <- fread(F_CELLS, colClasses = list(character = c("barcode", "removal_reason", "fail_reason_clean"))); cells4[is.na(removal_reason), removal_reason := ""]
basic <- readRDS(F_BASIC); MT <- basic[[2]]; stopifnot(length(MT) == 13, all(grepl("^MT-", MT)))
fsm <- fread(F_META1, select = c("V1", "ID", "name", "RNA_snn_res.3", "final_hit1_detailed", "final_hit1_group_detailed"), showProgress = FALSE)[ID %in% LIBS]
units <- fread(F_UNITS)[lib %in% LIBS]; nogrp <- readLines(F_NOGRP); dbl1 <- fread(F_DBL1)[lib %in% LIBS]
# carried label in the release vocabulary (carried_revised of helpers.R in 03_annotation/per_sample/)
source(file.path(NL, "scripts/helpers.R"), local = TRUE)

## ---- pass 1: per-library calls and fixed QC (no output per library yet) ----
P1 <- list(); MATS <- list()
for (L in LIBS) {
  x <- readRDS(file.path(TMPD, paste0(L, "_allbarcodes.rds"))); stopifnot(identical(x$lib, L), sum(x$cleaned@x) == man[lib == L, umi_clean])
  bc <- paste0(colnames(x$cleaned), "-1")
  if (L %in% NITERS_LIBS) {
    ed <- fread(file.path(TB, paste0(STUDY, "_02_ed_", L, "_niters.csv.gz")))
    call_set <- ed[called == TRUE, barcode]; fdr <- setNames(ed$FDR, ed$barcode); lim <- setNames(ed$Limited, ed$barcode)
    run_used <- paste0("emptyDrops(m, niters = ", NITERS_LIMITED, "), seed ", SEED_ED, " (re-call step 2)")
  } else {
    r <- rv[lib == L]; call_set <- r[called_cleaned_i == TRUE, barcode]; fdr <- setNames(r$FDR_cleaned_i, r$barcode); lim <- setNames(r$Limited_cleaned_i, r$barcode)
    run_used <- paste0("emptyDrops(m) DropletUtils defaults, seed ", SEED_ED, " (calls of the emptyDrops check, reused)")
  }
  stopifnot(all(call_set %in% bc), !anyDuplicated(call_set))
  if ("called_i_rule" %in% names(rv)) {   # other groups: the check's own rule column (cleaned_i, or cleaned_i_n30000 where Limited) must equal the call used here
    r2 <- rv[lib == L]; stopifnot(setequal(r2[called_i_rule == TRUE, barcode], call_set)); cat(L, ": call set equals the check's called_i_rule column\n") }
  c4 <- cells4[lib == L]; stopifnot(identical(c4$barcode, x$called))
  U <- union(x$called, call_set); j <- match(sub("-1$", "", U), colnames(x$cleaned)); stopifnot(!anyNA(j))
  m <- x$cleaned[, j, drop = FALSE]; colnames(m) <- U
  feat_name <- x$feat_name; stopifnot(length(feat_name) == nrow(m)); mtrow <- feat_name %in% MT; stopifnot(sum(mtrow) == 13)
  d <- data.table(study = STUDY, library = L, barcode = U, cell = paste0(L, "_", U), cr_called = U %in% x$called, recalled = U %in% call_set,
                  ed_FDR = unname(fdr[U]), ed_Limited = unname(lim[U]), ed_call_run = run_used,
                  umi_orig = unname(x$umi_orig[sub("-1$", "", U)]), umi_swapped = unname(x$umi_swapped[sub("-1$", "", U)]),
                  umi_clean = as.numeric(colSums(m)), genes_clean = as.integer(colSums(m > 0)), mt_clean = as.numeric(colSums(m[mtrow, , drop = FALSE])))
  stopifnot(all(d$umi_orig == d$umi_clean + d$umi_swapped))
  d[, frac_umi_removed := ifelse(umi_orig > 0, umi_swapped / umi_orig, NA_real_)]
  d[, pct_mt_clean := 100 * mt_clean / pmax(umi_clean, 1)]
  d[, pass_qc_clean := genes_clean >= 250 & umi_clean >= 500 & pct_mt_clean < 20]
  crit <- cbind(ifelse(d$genes_clean < 250, "genes<250", ""), ifelse(d$umi_clean < 500, "umi<500", ""), ifelse(d$pct_mt_clean >= 20, "mt>=20%", ""))
  d[, fail_reason_clean := apply(crit, 1, function(z) paste(z[z != ""], collapse = ";"))]
  ## the Cell Ranger-called rows must equal the swap step's cell table (independent code path)
  k <- match(c4$barcode, d$barcode)
  stopifnot(all(d$umi_clean[k] == c4$umi_clean), all(d$genes_clean[k] == c4$genes_clean), all(d$mt_clean[k] == c4$mt_clean),
            all(d$umi_orig[k] == c4$umi_orig), all(d$pass_qc_clean[k] == c4$pass_qc_clean), all(d$fail_reason_clean[k] == c4$fail_reason_clean))
  P1[[L]] <- d; MATS[[L]] <- list(m = m[, d$recalled, drop = FALSE], feat_id = rownames(x$cleaned), feat_name = feat_name)
  cat(sprintf("%s: Cell Ranger-called %d | re-called %d (%s) | both %d | re-called QC-pass %d\n", L, sum(d$cr_called), sum(d$recalled), run_used,
              sum(d$cr_called & d$recalled), sum(d$recalled & d$pass_qc_clean)))
  rm(x, m); invisible(gc())
}

## ---- library rule on the re-called set, then the pre-stated stop (before any library is processed) ----
D4 <- rbindlist(lapply(LIBS, function(L) { d <- P1[[L]]
  data.table(study = STUDY, library = L, recalled = sum(d$recalled), recalled_pass_qc_clean = sum(d$recalled & d$pass_qc_clean),
             decision4_pass_recalled = sum(d$recalled & d$pass_qc_clean) >= 100,
             p8_cr_called_pass_qc_clean = sum(d$cr_called & d$pass_qc_clean), p8_library_rule_pass = cells4[lib == L, library_rule_pass[1]],
             ed_call_run = d$ed_call_run[1]) }))
print(D4)
hit <- D4[library %in% PRESTATED_STOP & decision4_pass_recalled == TRUE, library]
stoptxt <- c(paste("Pre-stated stop for", STUDY, "at", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
             paste("libraries of the stop list in this group:", paste(intersect(PRESTATED_STOP, LIBS), collapse = ",")),
             paste("of these, passing the library rule on the re-called set:", if (length(hit)) paste(hit, collapse = ",") else "none"),
             if (length(hit)) "RESULT: STOP. No library of the group is processed. Handling the passing libraries needs a recorded decision of the author." else
               if (length(intersect(PRESTATED_STOP, LIBS))) "RESULT: not triggered; processing continues." else "RESULT: not applicable (no library of the stop list in this group); processing continues.")
cat(stoptxt, sep = "\n")
if (length(hit)) { safe_fwrite(D4, file.path(TB, paste0(STUDY, "_03_decision4_recalled.csv")))
  safe_lines(stoptxt, file.path(TB, paste0(STUDY, "_03_prestated_stop.txt"))); quit(save = "no", status = 5) }

## ---- pass 2: classes and re-called count matrices ----
colsum_csc <- function(x, ptr) { cs <- c(0, cumsum(x)); cs[ptr[-1] + 1] - cs[ptr[-length(ptr)] + 1] }
ALL <- list(); H5 <- list()
for (L in LIBS) {
  check_disk(L); d <- P1[[L]]; lp <- D4[library == L, decision4_pass_recalled]
  c4 <- cells4[lib == L]; k <- match(d$barcode, c4$barcode)
  d[, in_atlas := !is.na(k) & c4$in_atlas[k] %in% TRUE][, in_raw_object := !is.na(k) & c4$in_raw_object[k] %in% TRUE]
  d[, pass_qc_orig := c4$pass_qc_orig[k]][, p8_removal_reason := ifelse(is.na(k), "", c4$removal_reason[k])]
  d[, p8_n2 := in_atlas & p8_removal_reason == ""]
  d[, decision4_pass_recalled := lp]
  f <- fsm[ID == L]; jf <- match(d$cell, f$V1)
  d[, in_first_submission_metacomb := !is.na(jf)]
  d[, fs_label_raw := f$final_hit1_group_detailed[jf]][, fs_fine_raw := f$final_hit1_detailed[jf]][, fs_res3_cluster := as.character(f$`RNA_snn_res.3`[jf])]
  d[, carried_label := ifelse(p8_n2, carried_revised(fs_label_raw, fs_fine_raw), NA_character_)]
  stopifnot(all(d$in_first_submission_metacomb[d$in_atlas]))           # every atlas cell has a row in the cell metadata of the first atlas version
  u <- units[lib == L & set == "B"]; d[, setB_unit := u$unit_id[match(cell, u$cell)]]
  d[, no_group_raw_only := cell %in% nogrp]
  d[, fs_scDblFinder_class_step01 := dbl1$scDblFinder.class[match(cell, dbl1$cell)]]
  if (any(d$in_atlas & !is.na(d$fs_scDblFinder_class_step01))) stopifnot(all(d$fs_scDblFinder_class_step01[d$in_atlas & !is.na(d$fs_scDblFinder_class_step01)] == "singlet"))
  d[, class := fifelse(p8_n2 & recalled & lp, "A",
               fifelse(p8_n2 & !recalled, "B",
               fifelse(p8_n2 & recalled & !lp, "B2",
               fifelse(in_atlas & !p8_n2, "C",
               fifelse(cr_called & !in_atlas, "D",
               fifelse(!cr_called & recalled, "E", "?"))))))]
  if (L == "L071") {   # L071 (no cells in the first atlas version) is processed whole on its re-called capture
    nl71 <- fread(file.path(NL, "L071/cells_L071.csv.gz"), select = c("cell", "retained")); cat("input", file.path(NL, "L071/cells_L071.csv.gz"), "md5", md5(file.path(NL, "L071/cells_L071.csv.gz")), "\n")
    stopifnot(!any(d$in_atlas), all(nl71$cell %in% d$cell[d$cr_called]))
    d[, nl_L071_retained := cell %in% nl71[retained == TRUE, cell]]
    d[, class := fifelse(recalled, "F", "G")]
  } else d[, nl_L071_retained := NA]
  stopifnot(!any(d$class == "?"), all(d$cr_called | d$recalled))
  first_fate <- function(r) paste(c(if (r$in_raw_object) "in the first-submission raw object (not an atlas cell)" else "not in the first-submission raw object",
    if (!is.na(r$setB_unit)) paste0("set-B remove unit ", r$setB_unit, " (route (ii) only if restored by the curator)"),
    if (r$no_group_raw_only) "raw-object-only Luminal_No_Group cell (route (iii), not admitted)",
    if (isFALSE(r$pass_qc_orig)) "failed the fixed QC on the Cell Ranger counts",
    if (!is.na(r$fs_scDblFinder_class_step01)) paste0("first-submission step01 scDblFinder class ", r$fs_scDblFinder_class_step01)), collapse = "; ")
  dd <- which(d$class == "D"); d[, first_submission_fate := ""]
  if (length(dd)) d$first_submission_fate[dd] <- vapply(dd, function(q) first_fate(d[q]), "")
  ALL[[L]] <- d
  ## re-called count matrix (cleaned counts, re-called barcodes), 10x HDF5 v3 as 04_swappeddrops_study.R F1/F-A
  mm <- MATS[[L]]; cl <- mm$m; stopifnot(identical(colnames(cl), d$barcode[d$recalled]))
  h5 <- file.path(OD, "counts", paste0(L, "_swappeddrops_cleaned_recalled_feature_bc_matrix.h5"))
  if (file.exists(h5)) stop("refusing to overwrite ", h5)
  H5D <- file.path(TMPD, "h5_stage"); dir.create(H5D, showWarnings = FALSE); tmp <- file.path(H5D, paste0(L, ".tmp.h5")); if (file.exists(tmp)) file.remove(tmp)
  fh <- lt[ID == L, filtered_h5]; gen <- unique(as.character(h5read(fh, "matrix/features/genome"))); ra0 <- h5readAttributes(fh, "/"); h5closeAll()
  stopifnot(length(gen) == 1, identical(as.character(h5read(fh, "matrix/features/id")), mm$feat_id)); h5closeAll()
  chem <- fromJSON(h5read(lt[ID == L, molecule_info], "metrics_json"))$chemistry_description; h5closeAll()
  write10xCounts(tmp, cl, barcodes = colnames(cl), gene.id = mm$feat_id, gene.symbol = mm$feat_name, type = "HDF5", version = "3", genome = gen, overwrite = TRUE)
  fid <- H5Fopen(tmp); gid <- H5Gopen(fid, "/")
  for (nm in c("chemistry_description", "library_ids")) if (H5Aexists(gid, nm)) H5Adelete(gid, nm)
  h5writeAttribute(chem, h5obj = gid, name = "chemistry_description", variableLengthString = TRUE, asScalar = TRUE, encoding = "UTF-8")
  h5writeAttribute(L, h5obj = gid, name = "library_ids", variableLengthString = TRUE, asScalar = TRUE, encoding = "UTF-8")
  H5Gclose(gid); H5Fclose(fid)
  rb <- as.character(h5read(tmp, "matrix/barcodes")); rp <- as.numeric(h5read(tmp, "matrix/indptr")); rd <- as.numeric(h5read(tmp, "matrix/data"))
  rid <- as.character(h5read(tmp, "matrix/features/id")); ri <- as.integer(h5read(tmp, "matrix/indices")); ra <- h5readAttributes(tmp, "/"); h5closeAll()
  stopifnot(identical(rb, colnames(cl)), identical(rid, mm$feat_id), identical(ri, cl@i), all(rd == cl@x), all(colsum_csc(rd, rp) == d$umi_clean[d$recalled]),
            all(diff(rp) == d$genes_clean[d$recalled]), identical(as.character(ra$chemistry_description), chem), identical(as.character(ra$library_ids), L))
  stopifnot(system2("sync", tmp) == 0)
  H5[[L]] <- data.table(study = STUDY, library = L, recalled_barcodes = ncol(cl), umi_clean = sum(cl@x), h5 = h5, md5 = md5(tmp), bytes = file.size(tmp),
                        chemistry = chem, genome = gen, staged = tmp)
  cat(sprintf("%s: classes %s | staged h5 %s md5 %s\n", L, paste(names(table(d$class)), table(d$class), sep = "=", collapse = " "), tmp, md5(tmp)))
}
ALL <- rbindlist(ALL); H5 <- rbindlist(H5)
## group-level asserts against the swap step's selections
stopifnot(sum(ALL$cr_called) == nrow(cells4), sum(ALL$p8_n2) == nrow(cells4[in_atlas & removal_reason == ""]),
          sum(ALL$class == "C") == nrow(cells4[in_atlas & removal_reason != ""]))
## all asserts passed: move the staged h5 files (copy -> md5 equal -> remove the scratch copy) and write the tables
for (q in seq_len(nrow(H5))) { check_disk("h5 copy"); dst <- H5$h5[q]; part <- paste0(dst, ".tmp_", Sys.getpid())
  stopifnot(!file.exists(dst), file.copy(H5$staged[q], part), system2("sync", part) == 0, md5(part) == H5$md5[q], file.rename(part, dst), md5(dst) == H5$md5[q])
  file.remove(H5$staged[q]); cat("[write]", dst, "md5", md5(dst), "\n") }
H5[, staged := NULL]
safe_fwrite(D4, file.path(TB, paste0(STUDY, "_03_decision4_recalled.csv")))
safe_lines(stoptxt, file.path(TB, paste0(STUDY, "_03_prestated_stop.txt")))
safe_fwrite(ALL, file.path(TB, paste0(STUDY, "_03_barcodes_stage1.csv.gz")))
safe_fwrite(H5, file.path(TB, paste0(STUDY, "_03_recalled_count_matrices.csv")))
sm <- ALL[, .(cr_called = sum(cr_called), recalled = sum(recalled), A_route_i_candidates = sum(class == "A"), B_not_recalled = sum(class == "B"),
              B2_library_rule = sum(class == "B2"), C_p8_removed = sum(class == "C"), D_cr_non_atlas = sum(class == "D"), E_newly_called = sum(class == "E"),
              E_pass_qc_clean = sum(class == "E" & pass_qc_clean), F_L071_recalled = sum(class == "F"), G_L071_not_recalled = sum(class == "G")), by = library]
sm <- rbind(sm, cbind(library = "total", sm[, lapply(.SD, sum), .SDcols = -1]))
print(sm); safe_fwrite(sm, file.path(TB, paste0(STUDY, "_03_class_counts.csv")))
cat("ALL_ASSERTS_PASSED:", STUDY, "step 3\nend", format(Sys.time()), "\n"); print(sessionInfo())
