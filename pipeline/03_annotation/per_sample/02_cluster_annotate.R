# Code release header: clusters the singlets of one sample and labels each cluster with SingleR on cluster pseudobulks. In the run, it processed three samples: L071, L111 and L126.
# Input: the output of 01_qc_doublets.R and five SingleR annotation references (two paths are placeholders). Methods: Cell-type annotation.
# Output: <LIB>_step02_cluster_annotate.rds.
# Settings: sc.proc(pc = 30, res = 3), then get.psk and process_singleR per reference. The prostate reference gives the label, mapped to the atlas cell types.
# sc.proc, get.psk and process_singleR are in verbatim_functions.R, and the map and marker functions in helpers.R, both in this folder.
# get_lib and the run settings come from libraries.R and 00_env.R (not included).
# In 14 samples, a step not included re-clustered clusters whose top SingleR label was luminal, basal, club or hillock and relabeled them with SingleR.
# The released cell types also include the curator's cluster decisions (not included).
# 02_cluster_annotate.R <LIB> -- per-sample clustering, cluster-pseudobulk SingleR (5 references), hit1, 17-type map, markers.
#   sc.proc: NormalizeData, FindVariableFeatures(vst, 2000), ScaleData, RunPCA(30), RunUMAP(1:30), FindNeighbors(1:30), FindClusters(res 3)
#            [function-default seeds: RunPCA 42, RunUMAP 42, FindClusters 0]
#   get.psk(anno = "RNA_snn_res.3", sample = "name"): clusters with more than 5 cells, AggregateExpression sums, TMM, log2(CPM + 1)
#   process_singleR(fine.tune = TRUE) per reference -> cell-level pred_<ref> (label.final) and hit1_<ref> by seurat_clusters
#   prostate-reference hit1 and pred of a cluster without a pseudobulk -> "lowqual" -> final_hit1 = hit1_prostateref
#   mapped_label = map17_revised(final_hit1): group_detailed of the hit1 row in celltype_colorcode_update_17papers.xlsx (this folder), then
#   Adipocytes -> Endothelial cell, Neuroendocrine -> Luminal, Neutrophils -> Neutrophil, and Granulocytes and hit1 mast -> Mast cell
source("<PROJECT_DIR>/new_libraries_v1/scripts/00_env.R")
source(file.path(OUT, "scripts/libraries.R")); source(file.path(OUT, "scripts/verbatim_functions.R")); source(file.path(OUT, "scripts/helpers.R"))
suppressPackageStartupMessages({library(Seurat); library(Matrix); library(edgeR); library(SingleR)})
LIB <- commandArgs(TRUE)[1]; L <- get_lib(LIB)
env_report(); check_disk()
s1 <- readRDS(file.path(SCR, LIB, paste0(LIB, "_step01_qc_doublets.rds")))
stopifnot(s1$lib_pass, identical(s1$h5_md5, unname(tools::md5sum(L$h5))))
cat("singlets:", ncol(s1$cts.sub), " genes:", nrow(s1$cts.sub), " any MT- gene:", any(grepl("^MT-", rownames(s1$cts.sub))), "\n")

obj <- CreateSeuratObject(counts = s1$cts.sub)
obj$ID <- LIB; obj$sample_id <- L$sample_id; obj$name <- L$name  # metadata (name = sample_id_GleasonGrade, the pseudobulk sample key)
obj <- sc.proc(obj, pc = 30, res = 3)
cat("clusters (RNA_snn_res.3):", length(levels(obj$RNA_snn_res.3)), "\n"); print(table(obj$RNA_snn_res.3))
stopifnot(identical(as.character(obj$RNA_snn_res.3), as.character(obj$seurat_clusters)))

psk <- get.psk(obj = obj, anno = "RNA_snn_res.3", sample = "name")
cat("pseudobulk columns:", ncol(psk$dgelist), " genes:", nrow(psk$dgelist), "\n")
refpath <- c(prostateref = "<input: SingleR reference object prostate.new.ref.rds of the first-version annotation>",
             metref      = "<input: SingleR reference object prostate.met.ref.rds of the first-version annotation>",
             p1          = file.path(C, "p1.psk.ref.rds"),
             m2          = file.path(C, "m2.psk.ref.rds"),
             m2_group    = file.path(C, "m2.psk.annogroup.ref.rds"))
pred <- list()
for (r in names(refpath)) {
  ref <- readRDS(refpath[[r]]); cat("reference", r, "md5", tools::md5sum(refpath[[r]]), "\n")
  pred[[r]] <- process_singleR(psk, ref = ref, labels = ref)
}
for (r in names(pred)) {                                          # cell-level pred_<ref> and hit1_<ref> by cluster
  sm <- pred[[r]]$summary
  m <- setNames(sm$label.final, sm$anno); obj@meta.data[[paste0("pred_", r)]] <- unname(m[as.character(obj$seurat_clusters)])
  m <- setNames(sm$hit1, sm$anno);        obj@meta.data[[paste0("hit1_", r)]] <- unname(m[as.character(obj$seurat_clusters)])
}
obj$pred_prostateref[is.na(obj$pred_prostateref)] <- "lowqual"     # clusters of 5 cells or fewer have no pseudobulk
obj$hit1_prostateref[is.na(obj$hit1_prostateref)] <- "lowqual"
obj$final_hit1 <- obj$hit1_prostateref                            # the prostate reference gives the label
cm <- read_cellmatch()
obj$mapped_label <- map17_revised(obj$final_hit1, cm)             # atlas cell type (17-type map, see header)
obj$mapped_fine  <- fine_revised(obj$final_hit1, cm)              # finer type: detail_reduce of the same row, with the Granulocytes and Neutrophils edits
obj$mapped_label_original_route <- map17_original(obj$final_hit1, cm)   # recorded only: Neutrophils -> Monocyte, no mast or Granulocytes edit
print(table(obj$mapped_label, useNA = "ifany"))

# cluster-level SingleR summary (all 5 references)
cl_tab <- data.frame(cluster = levels(obj$RNA_snn_res.3), n_cells = as.integer(table(obj$RNA_snn_res.3)), stringsAsFactors = FALSE)
for (r in names(pred)) {
  sm <- pred[[r]]$summary; i <- match(cl_tab$cluster, sm$anno)
  cl_tab[[paste0(r, "_hit1")]] <- sm$hit1[i]; cl_tab[[paste0(r, "_hit1_score")]] <- round(as.numeric(sm$hit1.score[i]), 4)
  cl_tab[[paste0(r, "_hit2")]] <- sm$hit2[i]; cl_tab[[paste0(r, "_hit2_score")]] <- round(as.numeric(sm$hit2.score[i]), 4)
  cl_tab[[paste0(r, "_finetuned_label")]] <- sm$label[i]; cl_tab[[paste0(r, "_label_final")]] <- sm$label.final[i]
  cl_tab[[paste0(r, "_pruned_label")]] <- sm$label.prune[i]
}

# markers: canonical panel + audit panel (long table), top cluster markers
cnt <- LayerData(obj, assay = "RNA", layer = "counts"); dat <- LayerData(obj, assay = "RNA", layer = "data")
cl <- as.character(obj$RNA_snn_res.3)
pan_genes <- unique(unlist(PANEL)); audit <- read_audit_panel()
cat("panel genes:", length(pan_genes), "; absent from the matrix:", setdiff(pan_genes, rownames(cnt)), "\n")
cat("audit panel types:", length(audit), "; audit genes absent:", setdiff(unique(unlist(audit)), rownames(cnt)), "\n")
all_genes <- unique(c(pan_genes, unlist(audit)))
mstat <- cluster_marker_stats(cnt, dat, cl, all_genes)
mstat$panel_type <- sapply(mstat$gene, function(g) paste(c(names(PANEL)[sapply(PANEL, function(v) g %in% v)]), collapse = "|"))
mstat$panel_source <- ifelse(mstat$gene %in% names(PANEL_SOURCE), PANEL_SOURCE[mstat$gene], "")
mstat$audit_panel_type <- sapply(mstat$gene, function(g) paste(names(audit)[sapply(audit, function(v) g %in% v)], collapse = "|"))
Idents(obj) <- "RNA_snn_res.3"
mk <- FindAllMarkers(obj, only.pos = TRUE, verbose = FALSE)        # Seurat 5.1.0 defaults
cat("FindAllMarkers rows:", nrow(mk), "\n")
safe_save_rds(list(obj = obj, pred = pred, psk = psk, cl_tab = cl_tab, mstat = mstat, markers = mk),
              file.path(SCR, LIB, paste0(LIB, "_step02_cluster_annotate.rds")))
sessionInfo()
