# =============================================================================
# 01_QC_and_annotation.R
# =============================================================================
# Step 1 of the pipeline: Per-sample quality control, doublet detection,
# and pseudobulk-based cell type annotation for all 17 datasets.
#
# Cho et al. — A harmonized single-cell RNA-seq atlas of the human prostate
# Scientific Data (in submission)
#
# Workflow position:
#   [THIS FILE] → 02_data_merge_and_filter.py
#
# What this script does:
#   For each of the 163 samples across 17 studies:
#     1. Load 10x Genomics h5 count matrix (Cell Ranger v7.1)
#     2. Apply QC filters: genes >= 250, UMIs >= 500, MT% < 20
#     3. Detect and remove doublets (scDblFinder)
#     4. Cluster per sample at high resolution (res=3) for pseudobulk profiles
#     5. Annotate clusters with SingleR against prostate-specific references
#        (CellAnno resource expanded with Henry et al. normal prostate profiles)
#     6. Transfer cluster-level labels to individual cells
#     7. Merge all 17 studies into one Seurat object; export as h5ad
#        → also export: all_cells_prefilter_annotation.csv (cell annotation table used
#          by 02_data_merge_and_filter.py to filter the Python merged object)
#
# Key output files:
#   QC_summary_sample_level.csv — per-sample QC metrics (Supplementary Table 2)
#   all_cells_prefilter_annotation.csv       — cell-level annotation table
#   object.h5ad                     — merged Seurat object exported for Python
#
# Dependencies: Seurat (v5.1), SingleR, scDblFinder, edgeR, SeuratDisk
# =============================================================================

a <- readRDS("/mctp/share/users/hanbyul/projects/human_prostate_sc/cellranger_7.1/17datasets/Zenodo_file_submission/all_cells/all_cells_integrated.rds")
table(a$Paper_ID)
table(a$final_hit1_group_broad)
dim(a)
setwd("/mctp/share/users/hanbyul/projects/human_prostate_sc/cellranger_7.1/17datasets")

options(bitmapType = "cairo")
png("abc", type = "cairo"); dev.off()

library(openxlsx)
library(Seurat)
library(Matrix)
library(edgeR)
library(SingleR)
library(scales)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(scDblFinder)
library(tidyverse)
library(SeuratDisk)

# =============================================================================
# Color palettes
# =============================================================================
paperorder <- c(
  "P16_2018_CellRep_Henry", "P17_2022_JPathol_Joseph", "P19_2024_BioRxiv_TSC",
  "P14_2020_Elife_Crowley", "P15_2021_JPathol_Joseph", "P03_2020_MolCancer_Ma",
  "P01_2021_CellRep_Tuong", "P05_2023_NatCommun_Hirz", "P02_2022_NatCommun_Wong",
  "P10_2024_NatCommun_Graham", "P06_2022_MolCancer_Heidegger",
  "P11_2022_JExpMed_Masetti", "P21_2024_NPJPrecisOncol_Cheng",
  "P08_2020_CommunBiol_Dong", "P09_2022_EurUrol_Cheng",
  "P12_2021_CancerCell_Kfoury", "P07_2022_Science_Chan"
)

paper_col <- c(
  "P01"="#009CD2","P02"="#e6194B","P03"="#3cb44b","P05"="#ffe119","P06"="#4363d8",
  "P07"="#f58231","P08"="#911eb4","P09"="#46f0f0","P10"="#2b1700","P11"="#f032e6",
  "P12"="#bcf60c","P14"="#fabebe","P15"="#008080","P16"="#e6beff","P17"="#9A6324",
  "P19"="#FFFF33","P21"="#999999"
)

study_col <- setNames(
  c("#009CD2","#e6194B","#3cb44b","#ffe119","#4363d8","#f58231","#911eb4",
    "#46f0f0","#f032e6","#bcf60c","#fabebe","#008080","#e6beff","#9A6324",
    "#FFFF33","#999999","#2b1700"),
  c("P01_2021_CellRep_Tuong","P02_2022_NatCommun_Wong","P03_2020_MolCancer_Ma",
    "P05_2023_NatCommun_Hirz","P06_2022_MolCancer_Heidegger","P07_2022_Science_Chan",
    "P08_2020_CommunBiol_Dong","P09_2022_EurUrol_Cheng","P11_2022_JExpMed_Masetti",
    "P12_2021_CancerCell_Kfoury","P14_2020_Elife_Crowley","P15_2021_JPathol_Joseph",
    "P16_2018_CellRep_Henry","P17_2022_JPathol_Joseph","P19_2024_BioRxiv_TSC",
    "P21_2024_NPJPrecisOncol_Cheng","P10_2024_NatCommun_Graham")
)

GleasonGrade_col <- setNames(
  c("#71A45E","#87CEFA","#0047ab",
    "#FFEDA0","#FFEDA0","#FEB24C","#FEB24C","#FC4E2A","#E31A1C","#800026",
    "#9E9AC8","#6A51A3","#4A1486"),
  c("N","B","BPH","PCa","GG1","GG2","GG2_GG3","GG3","GG4","GG5",
    "mCRPC","mCRPC-NETD","mCRPC-teSCC")
)

tissue_col <- setNames(
  c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF","#999999"),
  c("Prostate","Bone","Liver","LymphNode","RetrocruralAbdomenMass",
    "PelvicSideWall","Mediastinal mass","Epidural","Brain")
)

# =============================================================================
# Helper functions
# =============================================================================

# Standard Seurat processing: normalize → HVG → scale → PCA → UMAP → cluster
sc.proc <- function(sc, n = 2000, pc = 30, res = 0.5,
                    normalize = TRUE, hvg = TRUE, pca = TRUE,
                    umap = TRUE, cluster = TRUE, assay = "RNA") {
  DefaultAssay(sc) <- assay
  if (normalize)  sc <- NormalizeData(sc, verbose = FALSE)
  if (hvg)        sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = n, verbose = FALSE)
  if (pca) {
    sc <- ScaleData(sc, verbose = FALSE)
    sc <- RunPCA(sc, npcs = pc, verbose = FALSE)
  }
  if (umap)    sc <- RunUMAP(sc, reduction = "pca", dims = 1:pc)
  if (cluster) {
    sc <- FindNeighbors(sc, dims = 1:pc)
    sc <- FindClusters(sc, resolution = res)
  }
  return(sc)
}

# DimPlot with per-group cell counts in legend
dim_manual <- function(sc, anno = "anno", col = NULL, labels = TRUE,
                       reduction = "umap", raster = NULL) {
  if (is.null(col)) {
    cls <- hue_pal()(length(unique(sc@meta.data[, anno])))
    names(cls) <- names(table(sc@meta.data[, anno]))
  } else {
    cls <- col[names(col) %in% names(table(sc@meta.data[, anno]))]
    sc@meta.data[, anno] <- factor(sc@meta.data[, anno], levels = names(cls))
  }
  DimPlot(sc, group.by = anno, label = labels, reduction = reduction, raster = raster) +
    scale_color_manual(
      labels = paste0(names(cls), " (", table(sc@meta.data[, anno])[names(cls)],
                      ", ", round(table(sc@meta.data[, anno])[names(cls)] /
                                    length(sc@meta.data[, anno]) * 100, 1), "%)"),
      values = cls
    ) + theme(aspect.ratio = 1)
}

# Pseudobulk aggregation (cluster × sample) — input for SingleR annotation
get.psk <- function(obj, anno, sample, keep.all.genes = FALSE, min.cell = 5, assay = "RNA") {
  obj$anno_sample <- paste(as.character(obj@meta.data[, anno]),
                           obj@meta.data[, sample], sep = ".")
  a <- table(obj$anno_sample)
  a <- a[a > min.cell]
  obj <- obj[, obj$anno_sample %in% names(a)]
  whole.psk <- Seurat::AggregateExpression(
    object = obj, assays = assay, return.seurat = FALSE, group.by = "anno_sample"
  )[[1]]
  colnames(whole.psk) <- names(a)
  whole.psk <- if (keep.all.genes) DGEList(counts = whole.psk) else
    DGEList(counts = whole.psk[rowSums(whole.psk) > 0, ])
  whole.psk <- calcNormFactors(whole.psk)
  logtpm <- log2(edgeR::cpm(whole.psk) + 1)
  s <- as.data.frame(t(sapply(colnames(whole.psk), function(x) strsplit(x, "[.]")[[1]])))
  colnames(s) <- c("anno", "sample")
  whole.psk$samples <- cbind.data.frame(whole.psk$samples, s)
  list(dgelist = whole.psk, logtpm = logtpm)
}

# =============================================================================
# Load metadata and gene lists
# =============================================================================
study.ids <- c(
  "P02_2022_NatCommun_Wong", "P03_2020_MolCancer_Ma", "P05_2023_NatCommun_Hirz",
  "P06_2022_MolCancer_Heidegger", "P07_2022_Science_Chan", "P08_2020_CommunBiol_Dong",
  "P09_2022_EurUrol_Cheng", "P11_2022_JExpMed_Masetti", "P12_2021_CancerCell_Kfoury",
  "P14_2020_Elife_Crowley", "P15_2021_JPathol_Joseph", "P16_2018_CellRep_Henry",
  "P17_2022_JPathol_Joseph", "P19_2024_BioRxiv_TSC", "P21_2024_NPJPrecisOncol_Cheng",
  "P10_2024_NatCommun_Graham", "P01_2021_CellRep_Tuong"
)

meta.clean  <- readRDS("meta.clean.update.rds")
# write.csv(meta.clean, "sample_metadata.csv", row.names = FALSE)
basic.info  <- readRDS("../basic.info.rds")
mt          <- basic.info[[2]]       # mitochondrial gene names (for PercentageFeatureSet)
ribo_genes  <- basic.info$ribo
hb_genes    <- basic.info$hb_genes

# Wong 2022: tumor/benign barcode assignments from HTO demultiplexing
wongbarcode <- readRDS("/mctp/share/users/yupingz/projects/aniket_public_prostate_analysis/cellranger7/wong.tumor.benign.barcodes.rds")

# =============================================================================
# Step 1: Per-sample QC and doublet detection
# =============================================================================
cat("\n===== Step 1: Per-sample QC and doublet detection =====\n")

for (xxx in 1:17) {
  study.name <- study.ids[xxx]
  k          <- which(meta.clean$Paper_ID == study.name)
  cat(sprintf("\nProcessing: %s (%d samples)\n", study.name, length(k)))

  obj <- lapply(k, function(i) {

    # Load 10x count matrix
    if (meta.clean$Paper_ID[i] == "P02_2022_NatCommun_Wong") {
      # Wong 2022: multi-sample HTO demultiplexed h5; subset to tumor or benign barcodes
      cts <- Read10X_h5(paste0(
        "/mctp/share/users/yupingz/cellranger/cellranger_7.1/counts/",
        meta.clean$countmatrixlocation[i], "/", meta.clean$patient_id[i],
        "/outs/per_sample_outs/", meta.clean$patient_id[i],
        "/count/sample_filtered_feature_bc_matrix.h5"
      ))$`Gene Expression`
      wb  <- wongbarcode[[meta.clean$patient_id[i]]]
      cts <- if (meta.clean$Group[i] == "PCa")
        cts[, row.names(wb)[wb$final.call == "T"]] else
        cts[, row.names(wb)[wb$final.call == "B"]]
    } else if (meta.clean$Paper_ID[i] == "P01_2021_CellRep_Tuong") {
      cts <- Read10X_h5(paste0(
        "/mctp/share/users/hanbyul/cellranger/cellranger_7.1/counts/",
        meta.clean$countmatrixlocation[i], "/", meta.clean$sample_id[i],
        "/outs/filtered_feature_bc_matrix.h5"
      ))
    } else {
      cts <- Read10X_h5(paste0(
        "/mctp/share/users/yupingz/cellranger/cellranger_7.1/counts/",
        meta.clean$countmatrixlocation[i], "/", meta.clean$sample_id[i],
        "/outs/filtered_feature_bc_matrix.h5"
      ))
    }

    # Prefix barcodes with sample ID to ensure uniqueness across all samples
    colnames(cts) <- paste(meta.clean$ID[i], colnames(cts), sep = "_")
    cts.ori       <- cts

    # Compute QC metrics
    d <- CreateSeuratObject(counts = cts)
    d <- NormalizeData(d)
    d[["percent.mt"]]   <- PercentageFeatureSet(d, features = mt)
    d[["percent.ribo"]] <- PercentageFeatureSet(d, features = ribo_genes)
    d[["percent.hb"]]   <- PercentageFeatureSet(d, features = hb_genes)

    if (sum(d$nFeature_RNA >= 250) <= 100) {
      cat(sprintf("  [LOW QUAL] %s: <100 cells pass gene filter\n", meta.clean$sample_id[i]))
      return(list(cts.ori = cts.ori, cts.filter = cts.ori, anno = cts.ori, cts.sub = cts.ori))
    }

    # Filter cells: >= 250 genes, >= 500 UMIs, < 20% MT content
    cts <- cts[, d$nFeature_RNA >= 250 & d$nCount_RNA >= 500 & d$percent.mt < 20]
    cts <- cts[!(rownames(cts) %in% mt), ]   # remove mitochondrial genes

    seu <- CreateSeuratObject(counts = cts)
    seu <- NormalizeData(seu)
    seu <- FindVariableFeatures(seu)
    seu <- ScaleData(seu)

    # Doublet detection (scDblFinder)
    sce     <- scDblFinder(as.SingleCellExperiment(seu))
    meta_dbl <- as.data.frame(sce@colData) %>% dplyr::select(starts_with("scDblFinder"))
    rownames(meta_dbl) <- rownames(sce@colData)
    seu <- AddMetaData(seu, meta_dbl %>% dplyr::select("scDblFinder.class"))

    doublet_anno <- meta_dbl %>% dplyr::select("scDblFinder.class")
    cts.sub      <- cts[, row.names(doublet_anno)[doublet_anno$scDblFinder.class == "singlet"]]

    list(cts.ori = cts.ori, cts.filter = cts, anno = doublet_anno, cts.sub = cts.sub)
  })

  saveRDS(obj, paste0("step01/", study.name, "_countmatrix_doublet.rds"))
}

# =============================================================================
# Step 2: Per-sample QC summary table (Supplementary Table 2)
# =============================================================================
cat("\n===== Step 2: QC summary =====\n")

qc_summary <- data.frame()

for (xxx in 1:17) {
  study.name <- study.ids[xxx]
  k   <- which(meta.clean$Paper_ID == study.name)
  obj <- readRDS(paste0("step01/", study.name, "_countmatrix_doublet.rds"))

  for (i in seq_along(obj)) {
    o <- obj[[i]]

    if (is.numeric(o$anno[1, 1])) {
      qc_summary <- rbind(qc_summary, data.frame(
        Study = study.name, Sample = meta.clean$sample_id[k[i]],
        ID = meta.clean$ID[k[i]], Exclude = meta.clean$Exclude[k[i]],
        Cells_before_QC = ncol(o$cts.ori),
        Median_Genes_before_QC = median(colSums(o$cts.ori > 0)),
        Median_UMIs_before_QC  = median(colSums(o$cts.ori)),
        Filtering_thresholds_used = "Genes>=250; UMIs>=500; MT%<20; remove MT genes",
        Cells_After_QC="LOW QUAL", Median_Genes_after_QC="LOW QUAL",
        Median_UMIs_after_QC="LOW QUAL", Cells_Doublet="LOW QUAL",
        Doublet_Rate_percent="LOW QUAL", Cells_Singlet="LOW QUAL",
        Median_Genes_only_singlet="LOW QUAL", Median_UMIs_only_singlet="LOW QUAL"
      ))
    } else {
      d <- CreateSeuratObject(counts = o$cts.ori)
      d <- NormalizeData(d)
      n_doublet    <- sum(o$anno == "doublet")
      doublet_rate <- round(n_doublet / nrow(o$anno) * 100, 2)

      qc_summary <- rbind(qc_summary, data.frame(
        Study = study.name, Sample = meta.clean$sample_id[k[i]],
        ID = meta.clean$ID[k[i]], Exclude = meta.clean$Exclude[k[i]],
        Cells_before_QC            = ncol(o$cts.ori),
        Median_Genes_before_QC     = median(colSums(o$cts.ori > 0)),
        Median_UMIs_before_QC      = median(colSums(o$cts.ori)),
        Percent_Mitochondrial_before_QC = round(median(PercentageFeatureSet(d, features = mt)), 2),
        Filtering_thresholds_used  = "Genes>=250; UMIs>=500; MT%<20; remove MT genes",
        Cells_After_QC             = ncol(o$cts.filter),
        Median_Genes_after_QC      = median(colSums(o$cts.filter > 0)),
        Median_UMIs_after_QC       = median(colSums(o$cts.filter)),
        Cells_Doublet              = n_doublet,
        Doublet_Rate_percent       = doublet_rate,
        Cells_Singlet              = ncol(o$cts.sub),
        Median_Genes_only_singlet  = median(colSums(o$cts.sub > 0)),
        Median_UMIs_only_singlet   = median(colSums(o$cts.sub))
      ))
    }
  }
}

qc_summary$GleasonScore <- meta.clean$GleasonScore[match(qc_summary$Sample, meta.clean$sample_id)]
qc_summary$GleasonGrade <- meta.clean$GleasonGrade[match(qc_summary$Sample, meta.clean$sample_id)]
qc_summary$patient_id   <- meta.clean$patient_id[match(qc_summary$Sample,  meta.clean$sample_id)]

saveRDS(qc_summary, "QC_summary_sample_level.rds")
write.csv(qc_summary, "QC_summary_sample_level.csv")
cat("Saved: QC_summary_sample_level.csv\n")

# =============================================================================
# Step 3: Per-sample clustering (res=3) for pseudobulk generation
# =============================================================================
cat("\n===== Step 3: Per-sample clustering and pseudobulk =====\n")

meta.clean <- readRDS("meta.clean.update.rds")
study.ids  <- names(study_col)

for (xxx in 1:17) {
  study.name  <- study.ids[xxx]
  include_ids <- meta.clean$ID[meta.clean$Paper_ID == study.name & meta.clean$Exclude == "Include "]
  all_ids     <- meta.clean$ID[meta.clean$Paper_ID == study.name]
  obj.ori     <- readRDS(paste0("step01/", study.name, "_countmatrix_doublet.rds"))
  names(obj.ori) <- all_ids
  obj.ori        <- obj.ori[names(obj.ori) %in% include_ids]

  obj          <- do.call(cbind, lapply(obj.ori, function(x) x$cts.sub))
  doublet_anno <- do.call(rbind, lapply(obj.ori, function(x) x$anno))
  obj          <- CreateSeuratObject(counts = obj)
  tmp          <- meta.clean[match(obj$orig.ident, meta.clean$ID), ]
  rownames(tmp)  <- rownames(obj@meta.data)
  tmp$doublet    <- as.character(doublet_anno[match(rownames(tmp), rownames(doublet_anno)), ])
  tmp$doublet[is.na(tmp$doublet)] <- "lowqual"
  obj <- AddMetaData(obj, tmp)
  obj$name <- paste0(obj$sample_id, "_", obj$GleasonGrade)

  obj.split <- SplitObject(obj, split.by = "name")
  obj.split <- lapply(obj.split, sc.proc, pc = 30, res = 3)
  saveRDS(obj.split, paste0("step02/", study.name, "_split_seurat_obj_singlet.rds"))
}

for (xxx in 1:17) {
  study.name <- study.ids[xxx]
  obj.split  <- readRDS(paste0("step02/", study.name, "_split_seurat_obj_singlet.rds"))
  psk        <- lapply(obj.split, get.psk, anno = "RNA_snn_res.3", sample = "name")
  saveRDS(psk, paste0("step02/", study.name, ".split.psk.rds"))
}

obj.split.psk.list <- lapply(1:17, function(i) readRDS(paste0("step02/", study.ids[i], ".split.psk.rds")))
names(obj.split.psk.list) <- sub("_GSE[0-9]+", "", study.ids)
saveRDS(obj.split.psk.list, "all.split.psk.list.rds")

# =============================================================================
# Step 4: SingleR annotation (pseudobulk, cluster-level)
# Reference: prostate-specific reference (CellAnno + Henry et al. prostate profiles)
# =============================================================================
cat("\n===== Step 4: SingleR annotation =====\n")

prostate.new.ref     <- readRDS("/mctp/share/users/yupingz/projects/aniket_public_prostate_analysis/cellranger7/prostate.new.ref.rds")
p1.psk.ref           <- readRDS("../p1.psk.ref.rds")
m2.psk.ref           <- readRDS("../m2.psk.ref.rds")
m2.psk.annogroup.ref <- readRDS("../m2.psk.annogroup.ref.rds")

process_singleR <- function(x, ref, labels) {
  p  <- SingleR(test = x[[2]], ref = ref[[1]], labels = labels[[2]], fine.tune = TRUE)
  s1 <- t(sapply(1:nrow(p$scores), function(i) names(sort(p$scores[i,], decreasing=TRUE)[1:2])))
  s2 <- t(sapply(1:nrow(p$scores), function(i) sort(p$scores[i,], decreasing=TRUE)[1:2]))
  s  <- as.data.frame(cbind(s1, s2))
  colnames(s) <- c("hit1","hit2","hit1.score","hit2.score")
  s$label       <- p$labels
  s$label.prune <- p$pruned.labels
  s$label.final <- ifelse(s$hit1 == s$label, s$hit1, NA)
  s <- cbind(x[[1]]$samples, s)
  list(result = p, summary = s)
}

obj.split.psk.list <- readRDS("all.split.psk.list.rds")

for (ref_name in c("prostate.new.ref","p1.psk.ref","m2.psk.ref","m2.psk.annogroup.ref")) {
  ref       <- get(ref_name)
  pred_list <- lapply(obj.split.psk.list, function(study) lapply(study, process_singleR, ref=ref, labels=ref))
  saveRDS(pred_list, paste0("obj.split.pred.list_", ref_name, ".rds"))
  cat(sprintf("Saved: %s predictions\n", ref_name))
}

# =============================================================================
# Step 5: Transfer cluster-level annotations to individual cells
# =============================================================================
cat("\n===== Step 5: Transfer annotations to cells =====\n")

obj.split.list <- readRDS("all.split.list.rds")

pred_prostateref <- readRDS("obj.split.pred.list_prostate.new.ref.rds")
pred_p1          <- readRDS("obj.split.pred.list_p1.psk.ref.rds")
pred_m2          <- readRDS("obj.split.pred.list_m2.psk.ref.rds")
pred_m2_group    <- readRDS("obj.split.pred.list_m2.psk.annogroup.ref.rds")

for (jj in 1:17) {
  obj.split <- obj.split.list[[jj]]
  for (i in seq_along(obj.split)) {
    for (pl in list(
      list(pred_prostateref[[jj]][[i]], "hit1_prostateref"),
      list(pred_p1[[jj]][[i]],          "hit1_p1"),
      list(pred_m2[[jj]][[i]],          "hit1_m2"),
      list(pred_m2_group[[jj]][[i]],    "hit1_m2_group")
    )) {
      m <- setNames(pl[[1]]$summary$hit1, pl[[1]]$summary$anno)
      obj.split[[i]][[pl[[2]]]] <-
        setNames(m[as.character(obj.split[[i]]$seurat_clusters)], colnames(obj.split[[i]]))
    }
  }
  obj.split.list[[jj]] <- obj.split
}

# =============================================================================
# Step 6: Merge all 17 studies → export as h5ad for Python pipeline
# The all_cells_prefilter_annotation.csv (cell annotation table) is exported here
# and used by 02_data_merge_and_filter.py to filter the Python merged object.
# =============================================================================
cat("\n===== Step 6: Merge and export =====\n")

obj.split.list.final <- obj.split.list

for (jj in 1:17) {
  obj.split <- obj.split.list[[jj]]
  cts <- do.call(cbind, lapply(obj.split, function(x) as.matrix(x@assays$RNA$counts)))
  tmp <- do.call(rbind, lapply(obj.split, function(x) x@meta.data))
  obj <- CreateSeuratObject(counts = cts)
  obj <- AddMetaData(obj, tmp)
  obj <- sc.proc(obj, pc = 30, res = 3)
  obj.split.list.final[[jj]] <- obj
}
names(obj.split.list.final) <- names(obj.split.list)
saveRDS(obj.split.list.final, "obj.split.list.cluster.rds")

pbmc.big <- merge(obj.split.list.final[[1]],
                  y = obj.split.list.final[2:17],
                  add.cell.ids = names(obj.split.list.final),
                  project = "prostate_atlas")

obj <- JoinLayers(pbmc.big)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$name)
obj <- sc.proc(obj, cluster = FALSE, pc = 30, res = 3)
saveRDS(obj, "obj.split.list.final.cluster.bypatient.rds")

# Export cell annotation table for Python (used in 02_data_merge_and_filter.py)
write.csv(obj@meta.data, "all_cells_prefilter_annotation.csv")
cat("Saved: all_cells_prefilter_annotation.csv\n")

# Export as h5ad
object <- obj
object[["RNA3"]] <- as(object[["RNA"]], Class = "Assay")
DefaultAssay(object) <- "RNA3"
object[["RNA"]] <- NULL
object <- RenameAssays(object, RNA3 = "RNA")
SaveH5Seurat(object, filename = "object.h5Seurat")
Convert("object.h5Seurat", dest = "h5ad")
cat("Saved: object.h5ad\n")
cat("→ Next: run 02_data_merge_and_filter.py\n")
