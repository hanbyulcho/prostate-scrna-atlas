# Code release header: a restore-mode variant of 04_doublets_cluster_annotate.R, for a sample of a barcode-swap study with restored cells.
# 04_doublets_cluster_annotate.R has its own restore mode. It writes to primary_b_v4/<STUDY>/<LIB>_restore_<RESTORE_TAG>/ and leaves the labels of restored cells blank.
# This variant takes those labels from RESTORE_LIST. Methods: Quality control, Cell-type annotation.
# Input: the inputs of 04_doublets_cluster_annotate.R and its cells_<LIB>.csv.gz of the same sample (run without RESTORE_LIST).
#   Environment variables: RESTORE_LIST and RESTORE_TAG (both required) and RECALL_INPUT_STATUS (required for every study except P09_Cheng).
#   Optional: RECALL_RUN_ROOT sets the run folder (its tables/ and outputs), and RECALL_TEST_OUT a test output folder.
# RESTORE_LIST names each restored cell (cell, library, setB_unit, label, answer_row_id).
# Each is a Cell Ranger-called cell outside the first atlas version, in a removed group (setB_unit) that the curator restored.
# Route ii keeps each listed cell that is re-called and passes cell QC, in a sample that passes the library rule. Each such cell takes the label RESTORE_LIST gives it.
# Every other step, seed and setting is that of 04_doublets_cluster_annotate.R. Outputs go to primary_b_v6/<STUDY>/<LIB>/ with its file names.
# Also written: cluster flags with the route-ii labels counted as carried labels (information only), and the route-v label and cluster changes.
# 04b_restored_cells.R STUDY LIB -- one sample (one GEM well). See 04_doublets_cluster_annotate.R for the steps.
source("<PROJECT_DIR>/recall_v1/scripts/00_env.R")
source(file.path(NL, "scripts/verbatim_functions.R")); source(file.path(NL, "scripts/helpers.R"))
suppressPackageStartupMessages({library(Seurat); library(Matrix); library(edgeR); library(SingleR); library(scDblFinder); library(SingleCellExperiment)
  library(dplyr); library(data.table)})
a <- commandArgs(TRUE); STUDY <- a[1]; LIB <- a[2]; stopifnot(STUDY %in% names(GROUPS))
script_header(); pkg_versions(); check_disk("start")
STATUS <- Sys.getenv("RECALL_INPUT_STATUS", ""); if (STUDY != "P09_Cheng") stopifnot(nzchar(STATUS)); cat("input_status:", STATUS, "\n")
lab <- function(x) { x <- as.data.table(x); x[, input_status := STATUS]; x }
LUMRERUN <- c(L154 = "D11_T_GG1", L056 = "SRR17518475_GG4")         # the barcode-swap samples among the 14 samples of the epithelial re-clustering
RESTORE_LIST <- Sys.getenv("RESTORE_LIST", ""); RESTORE_TAG <- Sys.getenv("RESTORE_TAG", "")
stopifnot(nzchar(RESTORE_LIST), file.exists(RESTORE_LIST), nzchar(RESTORE_TAG))   # restore mode required
cat("RESTORE mode:", RESTORE_LIST, "md5", md5(RESTORE_LIST), "tag", RESTORE_TAG, "\n")
stopifnot(md5(file.path(NL, "scripts/verbatim_functions.R")) == "35609a43c33765adf8123294db7914da")   # md5 of the copy that ran, which has a 2-line header above the released function lines
cat("verbatim_functions.R md5", md5(file.path(NL, "scripts/verbatim_functions.R")), "| helpers.R md5", md5(file.path(NL, "scripts/helpers.R")), "\n")
ROOT <- Sys.getenv("RECALL_RUN_ROOT", OUT); if (ROOT != OUT) { stopifnot(startsWith(ROOT, SCRB) || startsWith(ROOT, file.path(OUT, "runs"))); cat("RUN ROOT:", ROOT, "\n") }
OUTW <- Sys.getenv("RECALL_TEST_OUT", file.path(ROOT, "primary_b_v6")); if (OUTW != file.path(ROOT, "primary_b_v6")) { stopifnot(startsWith(OUTW, SCRB)); cat("TEST mode: outputs to", OUTW, "\n") }   # tests write to scratch only
TB <- file.path(ROOT, "tables"); OD <- file.path(OUTW, STUDY, LIB); dir.create(OD, recursive = TRUE, showWarnings = FALSE)   # output folder <STUDY>/<LIB>/
SD <- if (OUTW != file.path(ROOT, "primary_b_v6")) file.path(OUTW, "scratch_rds", STUDY, LIB) else file.path(SCRB, "primary_b_v6", basename(ROOT), STUDY, LIB); dir.create(SD, recursive = TRUE, showWarnings = FALSE)
for (f in list.files(SD, full.names = TRUE)) stop("scratch folder not empty, refusing to reuse: ", f)
for (f in list.files(OD, full.names = TRUE)) stop("output folder not empty, refusing to overwrite: ", f)

st <- fread(file.path(TB, paste0(STUDY, "_03_barcodes_stage1.csv.gz")), na.strings = c("NA", ""))[library == LIB]
st[is.na(p8_removal_reason), p8_removal_reason := ""]; st[is.na(first_submission_fate), first_submission_fate := ""]
n3 <- st$class == "D" & !is.na(st$setB_unit) & !st$recalled                                        # a set-B cell that is not re-called cannot be restored
st[n3, first_submission_fate := gsub("set-B remove unit (\\S+) \\(route \\(ii\\) only if restored by the curator\\)",
  "set-B remove unit \\1: not re-called by emptyDrops after swap removal, so it leaves and cannot be restored", first_submission_fate)]
stopifnot(all(grepl("cannot be restored", st$first_submission_fate[n3]))); cat("set-B cells not re-called (cannot be restored):", sum(n3), "\n")
H5T <- fread(file.path(TB, paste0(STUDY, "_03_recalled_count_matrices.csv")))[library == LIB]; stopifnot(nrow(H5T) == 1, md5(H5T$h5) == H5T$md5)
cat("inputs: stage 1", file.path(TB, paste0(STUDY, "_03_barcodes_stage1.csv.gz")), "md5", md5(file.path(TB, paste0(STUDY, "_03_barcodes_stage1.csv.gz"))), "\n")
D4 <- fread(file.path(TB, paste0(STUDY, "_03_decision4_recalled.csv")))[library == LIB]; lib_pass <- D4$decision4_pass_recalled
cat("count source:", H5T$h5, "md5", H5T$md5, "| library rule on the re-called set:", D4$recalled_pass_qc_clean, "->", lib_pass, "\n")
basic.info <- readRDS(file.path(C, "basic.info.rds")); mt <- basic.info[[2]]; ribo_genes <- basic.info$ribo; hb_genes <- basic.info$hb_genes
stopifnot(length(mt) == 13)

cts <- Read10X_h5(H5T$h5); if (is.list(cts)) cts <- cts$`Gene Expression`
stopifnot(all(grepl("-1$", colnames(cts))), !anyDuplicated(colnames(cts)), all(mt %in% rownames(cts)), all(cts@x == round(cts@x)))
raw_bc <- colnames(cts); colnames(cts) <- paste(LIB, raw_bc, sep = "_")
stopifnot(setequal(colnames(cts), st[recalled == TRUE, cell]))
d <- CreateSeuratObject(counts = cts)                                     # QC columns as in 01_qc_doublets.R
d <- NormalizeData(object = d, verbose = FALSE)
d[["percent.mt"]] <- PercentageFeatureSet(d, features = mt); d[["percent.ribo"]] <- PercentageFeatureSet(d, features = ribo_genes)
d[["percent.hb"]] <- PercentageFeatureSet(d, features = hb_genes)
k <- match(st$cell, colnames(d)); r <- !is.na(k)
st[, nCount_RNA_allgenes := NA_real_][, nFeature_RNA_allgenes := NA_real_][, pct_mt := NA_real_][, pct_ribo := NA_real_][, pct_hb := NA_real_]
st[r, nCount_RNA_allgenes := d$nCount_RNA[k[r]]][r, nFeature_RNA_allgenes := d$nFeature_RNA[k[r]]][r, pct_mt := d$percent.mt[k[r]]]
st[r, pct_ribo := d$percent.ribo[k[r]]][r, pct_hb := d$percent.hb[k[r]]]
## Seurat metrics == step-3 metrics (independent code paths)
stopifnot(all(st$nCount_RNA_allgenes[r] == st$umi_clean[r]), all(st$nFeature_RNA_allgenes[r] == st$genes_clean[r]),
          all(abs(st$pct_mt[r] - st$pct_mt_clean[r]) < 1e-9), all((st$nFeature_RNA_allgenes[r] >= 250 & st$nCount_RNA_allgenes[r] >= 500 & st$pct_mt[r] < 20) == st$pass_qc_clean[r]))
rm(d); invisible(gc())

out <- list(library = LIB, h5 = H5T$h5, h5_md5 = H5T$md5, lib_pass = lib_pass, seed = SEED_DBL)
st[, scDblFinder_score_seed0 := NA_real_][, scDblFinder_class_seed0 := NA_character_][, nCount_RNA_noMT := NA_real_][, nFeature_RNA_noMT := NA_real_]
if (lib_pass) {
  qcp <- st[recalled == TRUE & pass_qc_clean == TRUE, cell]
  cf <- cts[, qcp, drop = FALSE]; cf <- cf[!(rownames(cf) %in% mt), , drop = FALSE]           # remove the 13 MT genes
  seu <- CreateSeuratObject(counts = cf); seu <- NormalizeData(object = seu); seu <- FindVariableFeatures(object = seu); seu <- ScaleData(object = seu)
  sce <- as.SingleCellExperiment(seu); cat("SCE (whole re-called capture, fixed-QC pass, MT-free):", dim(sce), "\n")
  set.seed(SEED_DBL); sce <- scDblFinder(sce)                                                # defaults, seed 0
  sce2 <- as.SingleCellExperiment(seu); set.seed(SEED_DBL); sce2 <- scDblFinder(sce2)
  same_class <- identical(as.character(sce2$scDblFinder.class), as.character(sce$scDblFinder.class)); max_sd <- max(abs(sce2$scDblFinder.score - sce$scDblFinder.score))
  cat("scDblFinder re-run with the same seed: identical class:", same_class, "; max |score diff|:", max_sd, "\n"); rm(sce2)
  q <- match(st$cell, colnames(sce)); w <- !is.na(q)
  st[w, scDblFinder_score_seed0 := sce$scDblFinder.score[q[w]]][w, scDblFinder_class_seed0 := as.character(sce$scDblFinder.class)[q[w]]]
  st[w, nCount_RNA_noMT := Matrix::colSums(cf)[cell]][w, nFeature_RNA_noMT := Matrix::colSums(cf > 0)[cell]]
  out$dbl_rerun_identical <- same_class; out$dbl_rerun_max_score_diff <- max_sd; out$n_dbl_input <- ncol(sce)
  rm(seu, sce); invisible(gc())
} else { cat("LIBRARY FAILS the library rule on the re-called set: no doublet detection, no clustering.\n"); cf <- NULL }

## routes i and v; every condition is asserted on the columns (reverse check below)
st[, route := ""][, status := ""][, reason := ""]
st[class == "A", `:=`(route = "i", status = "kept", reason = "first-submission atlas cell re-called; label and first-submission doublet class carried")]
st[class == "B", `:=`(status = "removed", reason = "not re-called by emptyDrops after swap removal")]
st[class == "B2", `:=`(status = "removed", reason = "library_rule (on the re-called set)")]
st[class == "C", `:=`(status = "removed (swap step)", reason = paste0("removal reason of the swap step: ", p8_removal_reason))]
st[class == "D", `:=`(status = "first-submission fate kept", reason = paste0("Cell Ranger-called, not a first-submission atlas cell: ", first_submission_fate))]
st[class == "E" & !pass_qc_clean, `:=`(status = "not retained", reason = paste0("newly called; fails the fixed QC on the cleaned counts (", fail_reason_clean, ")"))]
st[class == "E" & pass_qc_clean & !decision4_pass_recalled, `:=`(status = "not retained", reason = "newly called; library fails the library rule on the re-called set")]
st[class == "E" & pass_qc_clean & decision4_pass_recalled & scDblFinder_class_seed0 %in% "doublet", `:=`(status = "not retained", reason = paste0("newly called; scDblFinder doublet (seed ", SEED_DBL, ")"))]
st[class == "E" & pass_qc_clean & decision4_pass_recalled & scDblFinder_class_seed0 %in% "singlet", `:=`(route = "v", status = "kept", reason = "newly called after swap removal: fixed QC, library rule, singlet")]
## L071: whole per-sample processing on the re-called capture; retained cells are route iv ('iv_L071')
st[class == "F" & !pass_qc_clean, `:=`(status = "not retained", reason = paste0("L071 re-called; fails the fixed QC on the cleaned counts (", fail_reason_clean, ")"))]
st[class == "F" & pass_qc_clean & !decision4_pass_recalled, `:=`(status = "not retained", reason = "L071 re-called; library fails the library rule on the re-called set")]
st[class == "F" & pass_qc_clean & decision4_pass_recalled & scDblFinder_class_seed0 %in% "doublet", `:=`(status = "not retained", reason = paste0("L071 re-called; scDblFinder doublet (seed ", SEED_DBL, ")"))]
st[class == "F" & pass_qc_clean & decision4_pass_recalled & scDblFinder_class_seed0 %in% "singlet", `:=`(route = "iv_L071", status = "kept", reason = "L071 re-called capture, whole new-library processing: fixed QC, library rule, singlet")]
st[class == "G" & nl_L071_retained %in% TRUE, `:=`(status = "removed", reason = "L071 route-(iv) cell of the new-library processing, not re-called by emptyDrops after swap removal")]
st[class == "G" & !(nl_L071_retained %in% TRUE), `:=`(status = "not retained", reason = "L071 Cell Ranger-called barcode not re-called; not retained by the new-library processing either")]
if (nzchar(RESTORE_LIST)) {   # (C) route ii: restored set-B cells, re-called, cell QC, library rule
  rld <- fread(RESTORE_LIST, colClasses = "character", na.strings = NULL); stopifnot(all(c("cell", "library", "setB_unit", "label", "answer_row_id") %in% names(rld)), !anyDuplicated(rld$cell))
  rl_all <- rld$cell; rld <- rld[library == LIB]; rl <- rld$cell   # this library's cells only (column library)
  stopifnot(all(startsWith(rl, paste0(LIB, "_"))), length(rl) > 0, all(nzchar(rld$label)))
  bad <- setdiff(rl, st$cell[st$class == "D" & !is.na(st$setB_unit)])
  if (length(bad)) stop("restore list holds cells of ", LIB, " that are not Cell Ranger-called set-B cells: ", paste(head(bad, 5), collapse = ","))
  mu <- match(rl, st$cell); if (!identical(st$setB_unit[mu], rld$setB_unit)) stop("restore list unit differs from the step-03 set-B unit")   # each listed cell must be in the same unit as in 03_recall_qc_classify.R
  cat("restore list:", length(rl_all), "cells in all;", length(rl), "of", LIB, "| units:", paste(unique(rld$setB_unit), collapse = ","), "\n")
  rr <- st$cell %in% rl & st$class == "D" & !is.na(st$setB_unit)
  rlab <- setNames(rld$label, rld$cell); rrow <- setNames(rld$answer_row_id, rld$cell)
  st[rr & recalled & pass_qc_clean & decision4_pass_recalled, `:=`(route = "ii", status = "kept",
     reason = paste0("set-B cell restored by the curator: re-called, fixed QC, library rule; unit ", setB_unit, " restored with label ", rlab[cell], " (row ", rrow[cell], " of the curator's answer table, not included)"))]
  st[rr & !(recalled & pass_qc_clean & decision4_pass_recalled), reason := paste0(reason, "; restored but not restorable here (re-called ", recalled, ", fixed QC ", pass_qc_clean, ")")]
  st[route == "ii", carried_label := unname(rlab[cell])]   # restored cells carry their unit's label
  stopifnot(all(st$setB_unit[st$route == "ii"] %in% rld$setB_unit), all(nzchar(st$carried_label[st$route == "ii"])))
  cat("restore mode: restored cells in this library", sum(rr), "| kept as route (ii)", sum(st$route == "ii"), "| not restorable here", sum(rr & st$route != "ii"), "\n")
  print(st[rr, .N, by = .(setB_unit, route, carried_label, recalled, pass_qc_clean)]) }
stopifnot(all(st$status != ""))
st[, doublet_class_applied := fifelse(route == "ii", "singlet (first-submission class carried: cell of the first-submission per-sample object; restored set-B cell)", fifelse(route == "i", "singlet (first-submission atlas cell; class carried, no score)",
                              fifelse(class %in% c("E", "F") & !is.na(scDblFinder_class_seed0), paste0(scDblFinder_class_seed0, " (scDblFinder seed 0, whole re-called capture)"), "")))]
st[, newly_called := !cr_called & recalled]
## reverse check: every barcode meeting the route-(v) conditions is route (v), and only those
cond_v <- st$newly_called & st$pass_qc_clean & st$decision4_pass_recalled & st$scDblFinder_class_seed0 %in% "singlet" & LIB != "L071"
cond_iv <- LIB == "L071" & st$recalled & st$pass_qc_clean & st$decision4_pass_recalled & st$scDblFinder_class_seed0 %in% "singlet"
stopifnot(identical(cond_iv, st$route == "iv_L071"))
stopifnot(identical(cond_v, st$route == "v"), all(!st$cr_called[st$route == "v"]))
cond_i <- st$p8_n2 & st$recalled & st$decision4_pass_recalled
stopifnot(identical(cond_i, st$route == "i"), all(st$pass_qc_clean[st$route == "i"]), all(st$in_atlas[st$route == "i"]))
kept <- st[status == "kept", cell]
## clustering input: (a) cells of the per-sample objects of the first atlas version, re-called, cell QC on the cleaned counts, library rule;
## (b) route v; (c) L071 route iv
stopifnot(!any(st$in_first_submission_metacomb & !is.na(st$fs_scDblFinder_class_step01) & st$fs_scDblFinder_class_step01 != "singlet"),   # the doublet record of the first atlas version agrees
          !any(st$in_first_submission_metacomb & st$pass_qc_orig %in% FALSE))
st[, clustered := recalled & pass_qc_clean & decision4_pass_recalled & (in_first_submission_metacomb %in% TRUE | route %in% c("v", "iv_L071"))]
stopifnot(all(st$clustered[st$status == "kept"]))
clin <- st[clustered == TRUE, cell]
cat(sprintf("clustering input: %d (kept %d: route (i) %d, (ii) %d, (v) %d, L071 %d; other singlets %d)\n", length(clin), length(kept), sum(st$route == "i"), sum(st$route == "ii"),
            sum(st$route == "v"), sum(st$route == "iv_L071"), sum(st$clustered & st$status != "kept")))
cat(sprintf("kept: route (i) %d | route (v) %d | L071 route (iv) %d | removed (not re-called) %d | newly called %d (QC-pass %d, doublets %d)\n", sum(st$route == "i"), sum(st$route == "v"), sum(st$route == "iv_L071"),
            sum(st$class == "B"), sum(st$newly_called), sum(st$newly_called & st$pass_qc_clean), sum(st$newly_called & st$scDblFinder_class_seed0 %in% "doublet")))

## clustering and SingleR on the clustered cells, as in 02_cluster_annotate.R
st[, cluster := NA_character_]
for (cc in c("hit1_prostateref", "pred_prostateref", "hit1_metref", "hit1_p1", "hit1_m2", "hit1_m2_group", "mapped_label", "mapped_fine")) st[, (cc) := NA_character_]
fsn <- fread(file.path(C17, "metacomb_17data_final.csv"), select = c("ID", "sample_id", "name"), showProgress = FALSE)[ID == LIB]
if (nrow(fsn)) { NAME <- unique(fsn$name); SID <- unique(fsn$sample_id); key_src <- "metacomb name"
} else { mcl <- fread(file.path(C17, "meta.clean.csv"), select = c("ID", "sample_id", "GleasonGrade"))[ID == LIB]
  NAME <- paste0(mcl$sample_id, "_", mcl$GleasonGrade); SID <- mcl$sample_id; key_src <- "meta.clean.csv sample_id_GleasonGrade (no metacomb rows)" }
if (lib_pass && length(kept) > 0) stopifnot(length(NAME) == 1, length(SID) == 1, !grepl("[.]", NAME))
cat("pseudobulk sample key of", LIB, ":", paste(NAME, collapse = ","), "(", key_src, ")\n")
tri <- NULL; lumtab <- NULL; alt <- NULL   # lumtab and alt are also defined for samples without clustering
if (lib_pass && length(kept) > 0) {
  obj <- CreateSeuratObject(counts = cf[, clin, drop = FALSE])                       # the clustered cells
  obj$ID <- LIB; obj$sample_id <- SID; obj$name <- NAME
  obj <- sc.proc(obj, pc = 30, res = 3)
  cat("clusters (RNA_snn_res.3):", length(levels(obj$RNA_snn_res.3)), "\n"); print(table(obj$RNA_snn_res.3))
  stopifnot(identical(as.character(obj$RNA_snn_res.3), as.character(obj$seurat_clusters)))
  psk <- get.psk(obj = obj, anno = "RNA_snn_res.3", sample = "name")
  refpath <- c(prostateref = "<input: SingleR reference object prostate.new.ref.rds of the first-version annotation>",
               metref      = "<input: SingleR reference object prostate.met.ref.rds of the first-version annotation>",
               p1 = file.path(C, "p1.psk.ref.rds"), m2 = file.path(C, "m2.psk.ref.rds"), m2_group = file.path(C, "m2.psk.annogroup.ref.rds"))
  pred <- list()
  for (rr in names(refpath)) { ref <- readRDS(refpath[[rr]]); cat("reference", rr, refpath[[rr]], "md5", md5(refpath[[rr]]), "\n")
    pred[[rr]] <- process_singleR(psk, ref = ref, labels = ref); rm(ref) }
  for (rr in names(pred)) { sm <- pred[[rr]]$summary
    m <- setNames(sm$label.final, sm$anno); obj@meta.data[[paste0("pred_", rr)]] <- unname(m[as.character(obj$seurat_clusters)])
    m <- setNames(sm$hit1, sm$anno); obj@meta.data[[paste0("hit1_", rr)]] <- unname(m[as.character(obj$seurat_clusters)]) }
  obj$pred_prostateref[is.na(obj$pred_prostateref)] <- "lowqual"; obj$hit1_prostateref[is.na(obj$hit1_prostateref)] <- "lowqual"
  obj$final_hit1 <- obj$hit1_prostateref
  cm <- read_cellmatch(); obj$mapped_label <- map17_revised(obj$final_hit1, cm); obj$mapped_fine <- fine_revised(obj$final_hit1, cm)
  obj$hit1_prostateref_res3 <- obj$hit1_prostateref; obj$mapped_label_res3 <- obj$mapped_label
  obj$triage_unit <- as.character(obj$RNA_snn_res.3); obj$lumcluster <- "rest"; lumtab <- NULL
  if (LIB %in% names(LUMRERUN)) {   # (B) epithelial re-clustering, in a child process (04l_lumrerun_child.R, not included)
    stopifnot(NAME == LUMRERUN[[LIB]])
    stopifnot(!grepl("leiden-fixed", Sys.getenv("LD_LIBRARY_PATH")))      # this process runs without the leiden-fixed library path; only the child gets it
    sub2 <- obj; Idents(sub2) <- "seurat_clusters"
    fin <- file.path(SD, paste0(LIB, "_lum_in.rds")); fout <- file.path(SD, paste0(LIB, "_lum_out.rds")); safe_rds(sub2, fin); rm(sub2)
    if (file.exists(fout)) stop("child output exists, refusing to reuse: ", fout)
    CH <- file.path(OUT, "scripts/04l_lumrerun_child.R"); cat("luminal re-run child", CH, "md5", md5(CH), "\n")
    rc <- system2("<R_HOME>/bin/Rscript", c("--no-environ", CH, fin, fout, NAME), stdout = file.path(SD, paste0(LIB, "_lum_child.log")), stderr = file.path(SD, paste0(LIB, "_lum_child.log")),
                  env = c(paste0("LD_LIBRARY_PATH=<LEIDEN_ENV>/lib:", Sys.getenv("LD_LIBRARY_PATH")), paste0("TMPDIR=", Sys.getenv("TMPDIR"))))
    cat(readLines(file.path(SD, paste0(LIB, "_lum_child.log"))), sep = "\n"); stopifnot(rc == 0, file.exists(fout))
    ch <- readRDS(fout); lu <- ch$lu; lpred <- ch$lpred; cat("leidenalg (child)", ch$leidenalg, "\n")
    if (ncol(lu) < 50) { lu$hit1_prostateref_rerun <- lu$hit1_prostateref; lu$pred_prostateref_rerun <- lu$pred_prostateref }
    lu$hit1_prostateref_rerun <- as.character(lu$hit1_prostateref_rerun); lu$hit1_prostateref_rerun[is.na(lu$hit1_prostateref_rerun)] <- "lowqual"
    jj <- match(colnames(lu), colnames(obj)); stopifnot(!anyNA(jj))
    obj$lumcluster[jj] <- as.character(lu$seurat_clusters)
    obj$hit1_prostateref_lumrerun <- as.character(obj$hit1_prostateref); obj$hit1_prostateref_lumrerun[jj] <- lu$hit1_prostateref_rerun
    obj$final_hit1 <- obj$hit1_prostateref_lumrerun                                           # the re-clustering hit1 becomes final_hit1
    obj$mapped_label <- map17_revised(obj$final_hit1, cm); obj$mapped_fine <- fine_revised(obj$final_hit1, cm)
    obj$triage_unit[jj] <- paste0("lum", as.character(lu$seurat_clusters))
    ## SingleR summary per re-clustered cluster (computed in the child with get.psk and process_singleR)
    if (ncol(lu) >= 50) {   # check: the same labels after the factor(levels = detailcelltype_col) step and NA -> lowqual
      chkl <- unname(setNames(lpred$prostateref$summary$hit1, lpred$prostateref$summary$anno)[as.character(lu$seurat_clusters)])
      dcn <- ch$dcn; chkl <- as.character(factor(chkl, levels = dcn[dcn %in% chkl])); chkl[is.na(chkl)] <- "lowqual"
      stopifnot(identical(chkl, unname(lu$hit1_prostateref_rerun)))
    } else lpred <- lapply(pred, function(p) { sm0 <- p$summary; list(summary = sm0[as.character(sm0$anno) %in% as.character(unique(lu$seurat_clusters)), , drop = FALSE]) })   # < 50 cells: no re-clustering, the res-3 calls stand
    lumtab <- data.table(cell = colnames(lu), lumcluster = as.character(lu$seurat_clusters), hit1_prostateref_rerun = lu$hit1_prostateref_rerun)
    cat("luminal re-run:", ncol(lu), "cells in", length(unique(lu$seurat_clusters)), "Leiden clusters\n"); print(table(lu$seurat_clusters, lu$hit1_prostateref_rerun))
    out$lum <- lu; out$lpred <- lpred
  }
  ulev <- unique(obj$triage_unit); ulev <- ulev[order(grepl("^lum", ulev), suppressWarnings(as.integer(sub("^lum", "", ulev))))]
  obj$triage_unit <- factor(obj$triage_unit, levels = ulev)
  cl_tab <- data.frame(cluster = levels(obj$triage_unit), n_cells = as.integer(table(obj$triage_unit)), stringsAsFactors = FALSE)
  cl_tab$unit_type <- ifelse(grepl("^lum", cl_tab$cluster), "luminal re-run Leiden cluster (as in the code of the first atlas version)", "res-3 cluster")
  cl_tab$res3_clusters <- sapply(cl_tab$cluster, function(u) comp_string(as.character(obj$RNA_snn_res.3[obj$triage_unit == u])))
  for (rr in names(pred)) { sm <- pred[[rr]]$summary; i <- match(cl_tab$cluster, sm$anno)
    if (!is.null(lumtab)) { sl <- lpred[[rr]]$summary; sl$anno <- paste0("lum", sl$anno); sm <- rbind(sm[, intersect(names(sm), names(sl))], sl[, intersect(names(sm), names(sl))]); i <- match(cl_tab$cluster, sm$anno) }
    cl_tab[[paste0(rr, "_hit1")]] <- sm$hit1[i]; cl_tab[[paste0(rr, "_hit1_score")]] <- round(as.numeric(sm$hit1.score[i]), 4)
    cl_tab[[paste0(rr, "_hit2")]] <- sm$hit2[i]; cl_tab[[paste0(rr, "_hit2_score")]] <- round(as.numeric(sm$hit2.score[i]), 4)
    cl_tab[[paste0(rr, "_finetuned_label")]] <- sm$label[i]; cl_tab[[paste0(rr, "_label_final")]] <- sm$label.final[i]; cl_tab[[paste0(rr, "_pruned_label")]] <- sm$label.prune[i] }
  md <- obj@meta.data; i <- match(st$cell, rownames(md)); w <- !is.na(i)
  stopifnot(setequal(rownames(md), clin))
  st[w, cluster := as.character(md$RNA_snn_res.3)[i[w]]]
  st[, triage_unit := NA_character_][, hit1_prostateref_res3 := NA_character_][, mapped_label_res3 := NA_character_][, lumcluster := NA_character_]
  st[w, triage_unit := as.character(md$triage_unit)[i[w]]][w, hit1_prostateref_res3 := as.character(md$hit1_prostateref_res3)[i[w]]]
  st[w, mapped_label_res3 := as.character(md$mapped_label_res3)[i[w]]][w, lumcluster := as.character(md$lumcluster)[i[w]]]
  for (cc in c("hit1_prostateref", "pred_prostateref", "hit1_metref", "hit1_p1", "hit1_m2", "hit1_m2_group", "mapped_label", "mapped_fine")) st[w, (cc) := as.character(md[[cc]])[i[w]]]
  ## markers as in 02_cluster_annotate.R: canonical panel + audit panel, FindAllMarkers defaults
  cnt <- LayerData(obj, assay = "RNA", layer = "counts"); dat <- LayerData(obj, assay = "RNA", layer = "data"); clv <- as.character(obj$triage_unit)
  pan_genes <- unique(unlist(PANEL)); audit <- read_audit_panel(); all_genes <- unique(c(pan_genes, unlist(audit)))
  mstat <- cluster_marker_stats(cnt, dat, clv, all_genes)
  mstat$panel_type <- sapply(mstat$gene, function(g) paste(c(names(PANEL)[sapply(PANEL, function(v) g %in% v)]), collapse = "|"))
  mstat$panel_source <- ifelse(mstat$gene %in% names(PANEL_SOURCE), PANEL_SOURCE[mstat$gene], "")
  mstat$audit_panel_type <- sapply(mstat$gene, function(g) paste(names(audit)[sapply(audit, function(v) g %in% v)], collapse = "|"))
  Idents(obj) <- "triage_unit"; mk <- FindAllMarkers(obj, only.pos = TRUE, verbose = FALSE); cat("FindAllMarkers rows:", nrow(mk), "\n")

  ## triage: carried = route i, new = routes v and iv_L071
  ret <- st[status == "kept"]; ret[, cluster := triage_unit]; cl_levels <- levels(obj$triage_unit)   # triage on the triage unit (= res-3 cluster except (B))
  tri <- data.frame(library = LIB, cluster = cl_levels, stringsAsFactors = FALSE)
  tri$n_cells <- sapply(cl_levels, function(k) sum(ret$cluster == k)); tri$n_route_i <- sapply(cl_levels, function(k) sum(ret$cluster == k & ret$route == "i"))
  tri$n_route_v <- sapply(cl_levels, function(k) sum(ret$cluster == k & ret$route == "v")); tri$n_route_iv_L071 <- sapply(cl_levels, function(k) sum(ret$cluster == k & ret$route == "iv_L071"))
  tri$n_route_ii <- sapply(cl_levels, function(k) sum(ret$cluster == k & ret$route == "ii"))
  stopifnot(all(tri$n_cells == tri$n_route_i + tri$n_route_v + tri$n_route_iv_L071 + tri$n_route_ii))
  tri$unit_type <- cl_tab$unit_type[match(cl_levels, cl_tab$cluster)]; tri$res3_clusters <- cl_tab$res3_clusters[match(cl_levels, cl_tab$cluster)]
  tri$n_new <- tri$n_route_v + tri$n_route_iv_L071; tri$pct_new <- ifelse(tri$n_cells > 0, round(100 * tri$n_new / tri$n_cells, 1), NA)
  tri$n_clustered <- sapply(cl_levels, function(k) sum(as.character(obj$triage_unit) == k)); tri$n_clustered_not_retained <- tri$n_clustered - tri$n_cells
  ul <- function(v, k) { x <- unique(as.character(v[as.character(obj$triage_unit) == k])); stopifnot(length(x) == 1); x }
  tri$automatic_label <- sapply(cl_levels, function(k) ul(obj$mapped_label, k)); tri$automatic_fine <- sapply(cl_levels, function(k) ul(obj$mapped_fine, k))
  maj <- lapply(cl_levels, function(k) { x <- ret$carried_label[ret$cluster == k & ret$route == "i"]
    if (!length(x)) return(list(lab = NA, share = NA, comp = "")); t <- table(x); top <- names(t)[t == max(t)]
    list(lab = paste(top, collapse = " / "), share = round(max(t) / length(x), 4), comp = comp_string(x), tie = length(top) > 1, set = top) })
  tri$carried_composition <- sapply(maj, `[[`, "comp"); tri$carried_majority <- sapply(maj, `[[`, "lab"); tri$carried_majority_share <- sapply(maj, `[[`, "share")
  tri$F0_no_first_submission_cells <- sum(ret$route == "i") == 0
  tri$F1_new_ge20pct <- tri$n_cells > 0 & tri$n_new >= 0.20 * tri$n_cells
  tri$F2_carried_conflict <- sapply(maj, function(m) !is.na(m$lab) && (isTRUE(m$tie) || m$share < 0.90))
  tri$F3_automatic_ne_carried <- sapply(seq_along(maj), function(q) !is.na(maj[[q]]$lab) && !(tri$automatic_label[q] %in% maj[[q]]$set))
  tri$F4_automatic_not_a_type <- tri$automatic_label %in% NOT_A_TYPE
  fl <- c("F0_no_first_submission_cells", "F1_new_ge20pct", "F2_carried_conflict", "F3_automatic_ne_carried", "F4_automatic_not_a_type")
  tri$flagged <- apply(tri[, fl], 1, any); tri$flag_reasons <- apply(tri[, fl], 1, function(v) paste(sub("_.*", "", fl[v]), collapse = "+"))
  ## curator review item: a flagged cluster with at least one new cell (the re-called clustering decides new cells only)
  tri$step_D_item <- tri$flagged & tri$n_new > 0
  tri$pseudobulk <- tri$n_clustered > 5   # the SingleR pseudobulk holds all clustered cells (get.psk, > 5 cells)
  st[, cluster_flagged := tri$flagged[match(triage_unit, tri$cluster)]]
  st[, candidate_label := NA_character_][, label_source := ""]
  st[route == "i", `:=`(candidate_label = carried_label, label_source = "carried (first submission)")]
  st[route %in% c("v", "iv_L071"), `:=`(candidate_label = mapped_label, label_source = ifelse(cluster_flagged, "automatic, pending review (cluster flagged)", "automatic (cluster not flagged)"))]
  st[status != "kept", label_source := "none (not kept)"]
  st[route == "ii", `:=`(candidate_label = carried_label, label_source = "carried: restored set-B unit label (the curator's answer, step 04 re-run)")]   # the unit's label from RESTORE_LIST
  ## information only: the same F2 / F3 rule with the route-ii labels in the carried composition (not used for any flag or label)
  maj2 <- lapply(cl_levels, function(k) { x <- ret$carried_label[ret$cluster == k & ret$route %in% c("i", "ii")]
    if (!length(x)) return(list(lab = NA, share = NA, comp = "")); t <- table(x); top <- names(t)[t == max(t)]
    list(lab = paste(top, collapse = " / "), share = round(max(t) / length(x), 4), comp = comp_string(x), tie = length(top) > 1, set = top) })
  alt <- tri[, c("library", "cluster", "n_cells", "n_route_i", "n_route_ii", "n_new", "automatic_label", "carried_composition", "carried_majority", "carried_majority_share",
                 "F0_no_first_submission_cells", "F1_new_ge20pct", "F2_carried_conflict", "F3_automatic_ne_carried", "F4_automatic_not_a_type", "flagged", "flag_reasons", "step_D_item")]
  alt$carried_composition_incl_ii <- sapply(maj2, `[[`, "comp"); alt$carried_majority_incl_ii <- sapply(maj2, `[[`, "lab"); alt$carried_majority_share_incl_ii <- sapply(maj2, `[[`, "share")
  alt$F2_incl_ii <- sapply(maj2, function(m) !is.na(m$lab) && (isTRUE(m$tie) || m$share < 0.90))
  alt$F3_incl_ii <- sapply(seq_along(maj2), function(q) !is.na(maj2[[q]]$lab) && !(tri$automatic_label[q] %in% maj2[[q]]$set))
  fl2 <- c("F0_no_first_submission_cells", "F1_new_ge20pct", "F2_incl_ii", "F3_incl_ii", "F4_automatic_not_a_type")
  alt$flagged_incl_ii <- apply(alt[, fl2], 1, any); alt$flag_reasons_incl_ii <- apply(alt[, fl2], 1, function(v) paste(sub("_.*", "", fl2[v]), collapse = "+"))
  alt$step_D_item_incl_ii <- alt$flagged_incl_ii & alt$n_new > 0
  alt$differs_from_primary <- alt$flagged_incl_ii != alt$flagged | alt$step_D_item_incl_ii != alt$step_D_item | alt$flag_reasons_incl_ii != alt$flag_reasons
  alt$note <- "information only: route-(ii) labels in the carried composition. The main triage (cluster_triage) counts them apart."
  cat("information only: clusters whose flags would differ with route-(ii) labels as carried composition:", sum(alt$differs_from_primary), "\n")

  ## review sheet (no study or author columns)
  ct <- cl_tab; q <- match(tri$cluster, ct$cluster)
  sheet <- tri[, c("library", "cluster", "unit_type", "res3_clusters", "step_D_item", "flagged", "flag_reasons", "n_cells", "n_route_i", "n_route_v", "n_route_iv_L071", "n_route_ii", "n_clustered", "n_clustered_not_retained", "pct_new")]
  sheet$automatic_label_17type <- tri$automatic_label; sheet$automatic_fine <- tri$automatic_fine
  sheet$prostateref_hit1 <- ct$prostateref_hit1[q]; sheet$prostateref_hit1_score <- ct$prostateref_hit1_score[q]
  sheet$prostateref_hit2 <- ct$prostateref_hit2[q]; sheet$prostateref_hit2_score <- ct$prostateref_hit2_score[q]
  sheet$prostateref_score_gap <- round(ct$prostateref_hit1_score[q] - ct$prostateref_hit2_score[q], 4)
  sheet$prostateref_finetuned_label <- ct$prostateref_finetuned_label[q]; sheet$prostateref_label_final <- ct$prostateref_label_final[q]
  for (rr in c("metref", "p1", "m2", "m2_group")) { sheet[[paste0(rr, "_hit1")]] <- ct[[paste0(rr, "_hit1")]][q]
    sheet[[paste0(rr, "_hit1_score")]] <- ct[[paste0(rr, "_hit1_score")]][q]; sheet[[paste0(rr, "_label_final")]] <- ct[[paste0(rr, "_label_final")]][q] }
  sheet$carried_composition <- tri$carried_composition; sheet$carried_majority <- tri$carried_majority; sheet$carried_majority_share <- tri$carried_majority_share
  med <- function(v, k, rt = c("i", "v", "iv_L071", "ii")) { x <- ret[[v]][ret$cluster == k & ret$route %in% rt]; if (length(x)) round(median(x), 1) else NA }
  sheet$median_nCount_allgenes <- sapply(sheet$cluster, function(k) med("nCount_RNA_allgenes", k)); sheet$median_nFeature_allgenes <- sapply(sheet$cluster, function(k) med("nFeature_RNA_allgenes", k))
  sheet$median_pct_mt <- sapply(sheet$cluster, function(k) med("pct_mt", k))
  sheet$median_nCount_allgenes_new <- sapply(sheet$cluster, function(k) med("nCount_RNA_allgenes", k, c("v", "iv_L071")))
  sheet$median_nFeature_allgenes_new <- sapply(sheet$cluster, function(k) med("nFeature_RNA_allgenes", k, c("v", "iv_L071")))
  sheet$median_scDblFinder_score_new <- sapply(sheet$cluster, function(k) { x <- ret$scDblFinder_score_seed0[ret$cluster == k & ret$route %in% c("v", "iv_L071")]; if (length(x)) round(median(x), 3) else NA })
  for (ty in names(PANEL)) sheet[[paste0("markers_", gsub(" ", "_", ty))]] <- sapply(sheet$cluster, function(k) fmt_type(mstat[mstat$cluster == k, ], PANEL[[ty]]))
  top_str <- function(k, min_pct1 = 0) { x <- mk[mk$cluster == k & mk$p_val_adj < 0.05 & mk$pct.1 >= min_pct1, ]
    if (!nrow(x)) return(""); x <- head(x[order(-x$avg_log2FC), ], 10)
    paste(sprintf("%s (%.2f; %.0f%%/%.0f%%)", x$gene, x$avg_log2FC, 100 * x$pct.1, 100 * x$pct.2), collapse = "; ") }
  sheet$top10_markers_by_log2FC <- sapply(sheet$cluster, top_str); sheet$top10_markers_pct1_ge25 <- sapply(sheet$cluster, top_str, min_pct1 = 0.25)
  sheet$reviewer_label <- ""; sheet$reviewer_note <- ""
  sheet <- sheet[order(!sheet$step_D_item, !sheet$flagged, match(sheet$cluster, cl_levels)), ]
  bg <- sheet[1, ]; bg[] <- NA; bg$library <- LIB; bg$cluster <- "ALL (library background, all kept cells)"; bg$n_cells <- nrow(ret)
  for (ty in names(PANEL)) bg[[paste0("markers_", gsub(" ", "_", ty))]] <- fmt_type(mstat[mstat$cluster == "ALL_CELLS", ], PANEL[[ty]])
  bg$reviewer_label <- ""; bg$reviewer_note <- ""; sheet_out <- rbind(sheet, bg)
  banned <- grep("author|Paper_ID|study|source_label|original_label", colnames(sheet_out), ignore.case = TRUE, value = TRUE)
  if (length(banned)) stop("author-blind check failed: ", paste(banned, collapse = ", "))
  out$obj <- obj; out$pred <- pred; out$psk <- psk; out$cl_tab <- cl_tab; out$mstat <- mstat; out$markers <- mk
} else { st[, cluster_flagged := NA][, candidate_label := NA_character_][, label_source := ifelse(status == "kept", "", "none (not kept)")]
  st[, triage_unit := NA_character_][, lumcluster := NA_character_][, hit1_prostateref_res3 := NA_character_][, mapped_label_res3 := NA_character_] }

## information-only comparison: route-i cells, carried label vs the automatic label of their re-called cluster
cmpd <- st[route == "i", .(cell, library, route, cluster, carried_label_first_submission = fs_label_raw, carried_fine_first_submission = fs_fine_raw,
                           carried_label, fs_res3_cluster, automatic_label = mapped_label, automatic_fine = mapped_fine, hit1_prostateref)]
cmpd[, same_label := carried_label == automatic_label]
cm <- read_cellmatch(); broad <- c(setNames(cm$group_broad, cm$group_detailed), "Mast cell" = "Immune Cells", "Neutrophil" = "Immune Cells"); broad <- broad[!duplicated(names(broad))]
cmpd[, carried_compartment := unname(broad[carried_label])][, automatic_compartment := unname(broad[automatic_label])][, same_compartment := carried_compartment == automatic_compartment]

## per-cell table
cols <- c("cell", "barcode", "library", "route", "status", "reason", "cr_called", "recalled", "newly_called", "ed_FDR", "ed_Limited", "ed_call_run",
  "umi_orig", "umi_swapped", "frac_umi_removed", "umi_clean", "genes_clean", "mt_clean", "pct_mt_clean", "pass_qc_clean", "fail_reason_clean",
  "pct_ribo", "pct_hb", "decision4_pass_recalled", "nCount_RNA_noMT", "nFeature_RNA_noMT", "scDblFinder_score_seed0", "scDblFinder_class_seed0", "doublet_class_applied",
  "in_atlas", "in_raw_object", "pass_qc_orig", "p8_removal_reason", "p8_n2", "setB_unit", "no_group_raw_only", "fs_scDblFinder_class_step01", "first_submission_fate", "nl_L071_retained",
  "fs_label_raw", "fs_fine_raw", "carried_label", "fs_res3_cluster", "class",
  "cluster", "triage_unit", "lumcluster", "hit1_prostateref_res3", "mapped_label_res3", "hit1_prostateref", "pred_prostateref", "hit1_metref", "hit1_p1", "hit1_m2", "hit1_m2_group", "mapped_label", "mapped_fine",
  "cluster_flagged", "candidate_label", "label_source", "clustered")
o <- st[, ..cols]; o[, pct_mt_clean := round(pct_mt_clean, 4)][, pct_ribo := round(pct_ribo, 4)][, pct_hb := round(pct_hb, 4)][, frac_umi_removed := round(frac_umi_removed, 6)]
o[, scDblFinder_score_seed0 := signif(scDblFinder_score_seed0, 6)]
setnames(o, c("fs_label_raw", "fs_fine_raw", "fs_res3_cluster"), c("carried_label_first_submission", "carried_fine_first_submission", "first_submission_res3_cluster"))
safe_fwrite(lab(o), file.path(OD, paste0("cells_", LIB, ".csv.gz")))
if (!is.null(tri)) {
  safe_fwrite(lab(tri), file.path(OD, paste0("cluster_triage_", LIB, ".csv")))
  safe_fwrite(lab(sheet_out), file.path(OD, paste0("cluster_review_", LIB, ".csv")))
  safe_fwrite(lab(alt), file.path(OD, paste0("cluster_triage_route_ii_as_carried_", LIB, ".csv")))   # information only
  ms <- as.data.table(mstat); ms[, library := LIB]
  safe_fwrite(lab(ms[, .(library, cluster, gene, panel_type, panel_source, audit_panel_type, in_matrix, pct_expressing, mean_lognorm)]), file.path(OD, paste0("cluster_marker_panel_", LIB, ".csv.gz")))
  tm <- as.data.table(mk)[p_val_adj < 0.05][order(as.integer(as.character(cluster)), -avg_log2FC)]; tm[, library := LIB]
  safe_fwrite(lab(tm[, .(library, cluster, gene, avg_log2FC, pct.1, pct.2, p_val, p_val_adj)]), file.path(OD, paste0("cluster_markers_padj05_", LIB, ".csv.gz")))
}
safe_fwrite(lab(cmpd), file.path(OD, paste0("comparison_route_i_", LIB, ".csv.gz")))
S <- data.table(library = LIB, count_source = H5T$h5, count_source_md5 = H5T$md5, barcodes_in_table = nrow(st), cr_called = sum(st$cr_called), recalled = sum(st$recalled),
  newly_called = sum(st$newly_called), newly_called_pass_qc = sum(st$newly_called & st$pass_qc_clean), decision4_recalled_pass_qc = D4$recalled_pass_qc_clean, decision4_pass = lib_pass,
  dbl_input_cells = if (lib_pass) out$n_dbl_input else NA, clustered_cells = sum(st$clustered), clustered_not_retained = sum(st$clustered & st$status != "kept"), dbl_doublets_all = sum(st$scDblFinder_class_seed0 %in% "doublet"),
  dbl_doublets_newly_called = sum(st$newly_called & st$scDblFinder_class_seed0 %in% "doublet"), dbl_doublets_route_i_not_applied = sum(st$class == "A" & st$scDblFinder_class_seed0 %in% "doublet"),
  dbl_rerun_identical = if (lib_pass) out$dbl_rerun_identical else NA,
  route_i = sum(st$route == "i"), route_ii = sum(st$route == "ii"), route_v = sum(st$route == "v"), luminal_rerun = !is.null(lumtab), route_iv_L071 = sum(st$route == "iv_L071"), L071_route_iv_not_recalled = sum(st$class == "G" & st$nl_L071_retained %in% TRUE), kept = sum(st$status == "kept"), removed_not_recalled = sum(st$class == "B"), removed_library_rule = sum(st$class == "B2"),
  p8_removed = sum(st$class == "C"), cr_non_atlas = sum(st$class == "D"), cr_non_atlas_recalled = sum(st$class == "D" & st$recalled),
  clusters = if (is.null(tri)) 0L else nrow(tri), clusters_le5 = if (is.null(tri)) 0L else sum(!tri$pseudobulk),
  clusters_flagged = if (is.null(tri)) 0L else sum(tri$flagged), clusters_with_route_v = if (is.null(tri)) 0L else sum(tri$n_route_v > 0),
  step_D_clusters = if (is.null(tri)) 0L else sum(tri$step_D_item), route_v_in_step_D_clusters = sum(st$route %in% c("v", "iv_L071") & st$cluster_flagged %in% TRUE),
  route_v_automatic_unflagged = sum(st$route %in% c("v", "iv_L071") & st$cluster_flagged %in% FALSE),
  comparison_same_label = sum(cmpd$same_label, na.rm = TRUE), comparison_n = nrow(cmpd), comparison_same_compartment = sum(cmpd$same_compartment, na.rm = TRUE))
if (!is.null(lumtab)) {   # (B) comparison of the re-clustering label with the res-3 label, per kept cell
  lc <- st[status == "kept", .(cell, library, route, res3_cluster = cluster, hit1_prostateref_res3, mapped_label_res3, lumcluster, triage_unit,
                               hit1_prostateref_lumrerun = hit1_prostateref, mapped_label_lumrerun = mapped_label)]
  lc[, same_mapped_label := mapped_label_res3 == mapped_label_lumrerun]
  safe_fwrite(lab(lc), file.path(OD, paste0("lumrerun_vs_prostateref_cells_", LIB, ".csv.gz")))
  lsum <- lc[route %in% c("v", "iv_L071"), .(n_route_v = .N, labels_prostateref_res3 = comp_string(mapped_label_res3), labels_lumrerun = comp_string(mapped_label_lumrerun),
                                              n_same = sum(same_mapped_label)), by = .(library, res3_cluster)]
  safe_fwrite(lab(lsum), file.path(OD, paste0("lumrerun_vs_prostateref_route_v_clusters_", LIB, ".csv")))
}
if (nzchar(RESTORE_LIST)) {   # (C) route-(v) label changes against the run without restorations
  base <- fread(file.path(ROOT, "primary_b_v4", STUDY, LIB, paste0("cells_", LIB, ".csv.gz")), select = c("cell", "route", "cluster", "candidate_label", "label_source"), colClasses = list(character = "cluster"))[route == "v"]   # the run without restorations
  now <- st[route == "v", .(cell, cluster_restore = cluster, candidate_label_restore = candidate_label, label_source_restore = label_source)]
  dv <- merge(base, now, by = "cell", all = TRUE); dv[, label_changed := candidate_label != candidate_label_restore][, source_changed := label_source != label_source_restore]
  dv[, cluster_changed := cluster != cluster_restore]; stopifnot(!anyNA(dv$cluster_changed)); cat("route-(v) cells whose cluster differs from the run without restorations:", sum(dv$cluster_changed), "\n")   # route-v cells whose cluster changed
  safe_fwrite(lab(dv), file.path(OD, paste0("route_v_label_changes_", LIB, "_restore_", RESTORE_TAG, ".csv")))
  cat("restore mode: route-(v) cells", nrow(dv), "| label changed", sum(dv$label_changed, na.rm = TRUE), "| label_source changed", sum(dv$source_changed, na.rm = TRUE), "\n")
}
safe_fwrite(lab(S), file.path(OD, paste0("summary_", LIB, ".csv")))
safe_rds(out, file.path(SD, paste0(LIB, "_step04.rds")))
print(t(S)); if (!is.null(tri)) print(tri[tri$flagged, c("cluster", "n_cells", "n_route_i", "n_route_v", "automatic_label", "carried_majority", "carried_majority_share", "flag_reasons", "step_D_item")])
cat("ALL_ASSERTS_PASSED:", STUDY, LIB, "step 4\nend", format(Sys.time()), "\n"); print(sessionInfo())
