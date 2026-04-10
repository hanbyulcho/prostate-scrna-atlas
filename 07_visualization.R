# =============================================================================
# 07_visualization.R
# =============================================================================
# Step 7 of the pipeline: all publication figures and supplementary panels.
#
# Cho et al. — A harmonized single-cell RNA-seq atlas of the human prostate
# Scientific Data (in submission)
#
# Workflow position:
#   06_h5ad_to_seurat.R → [THIS FILE]
#
# What this script produces:
#
#   Section 1 — QC figures (Supplementary Figure 1a–d)
#     Source data: QC_metrics_raw_prefilter.csv, qc_summary_forsupplementary.csv
#     1a. Violin plots: nCount, nFeature, pct_mt across 17 studies (pre-filter)
#     1b. Bar chart: cells before QC / after QC / singlets / after annotation QC
#     1c. Stacked bar: cells removed at each filtering step
#     1d. Doublet rate per sample (boxplot + jitter)
#
#   Section 2 — All-cells scANVI UMAP (Figure 1, Supplementary Figure 1b–c)
#     Source: merged_a17_hit1_scib_50pc2000g.rds
#     DimPlots colored by: cell type, study, Group, GleasonGrade
#
#   Section 3 — Luminal scANVI UMAP (Figure 2)
#     Source: seu_lum_v1_noNoGroup_scib_30pc3000g.rds
#     DimPlots colored by: etsgroup, Group, GleasonGrade, study
#
#   Section 4 — scANVI feature plots for luminal cells (Supplementary Figure 4d)
#     Source: seu_lum_v1_noNoGroup_scib_30pc3000g.rds
#     Genes: AR, ERG, ETV1, ETV4, SPINK1 + AR-TCGA score + NE-Beltran score
#
#   Section 5 — Pseudobulk bubble (dot) plots (Figure 4c–d)
#     Source: pseudobulk Seurat RDS (correctmetapca10_new/seurat_all_allgroups_pca10_0320.rds)
#     Groups: primary PCa (Benign / nonETS / ERG / ETV); mCRPC (AR / NEPC / DN)
#     Marker genes selected by FindAllMarkers + known biology
#
# Dependencies: Seurat (v5), ggplot2, dplyr, tidyr, patchwork, gridExtra, ggrastr
# =============================================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# =============================================================================
# Shared color palettes (consistent across all sections)
# =============================================================================

paperorder <- c(
  "P07_2022_Science_Chan",        "P12_2021_CancerCell_Kfoury",
  "P09_2022_EurUrol_Cheng",       "P08_2020_CommunBiol_Dong",
  "P21_2024_NPJPrecisOncol_Cheng","P11_2022_JExpMed_Masetti",
  "P06_2022_MolCancer_Heidegger", "P10_2024_NatCommun_Graham",
  "P02_2022_NatCommun_Wong",      "P05_2023_NatCommun_Hirz",
  "P01_2021_CellRep_Tuong",       "P03_2020_MolCancer_Ma",
  "P15_2021_JPathol_Joseph",      "P14_2020_Elife_Crowley",
  "P19_2024_BioRxiv_TSC",         "P17_2022_JPathol_Joseph",
  "P16_2018_CellRep_Henry"
)

study_labels <- c(
  "P16_2018_CellRep_Henry"        = "Henry2018",
  "P17_2022_JPathol_Joseph"       = "Joseph2022",
  "P19_2024_BioRxiv_TSC"          = "TSC2024",
  "P14_2020_Elife_Crowley"        = "Crowley2020",
  "P15_2021_JPathol_Joseph"       = "Joseph2021",
  "P03_2020_MolCancer_Ma"         = "Ma2020",
  "P01_2021_CellRep_Tuong"        = "Tuong2021",
  "P05_2023_NatCommun_Hirz"       = "Hirz2023",
  "P02_2022_NatCommun_Wong"       = "Wong2022",
  "P10_2024_NatCommun_Graham"     = "Graham2024",
  "P06_2022_MolCancer_Heidegger"  = "Heidegger2022",
  "P11_2022_JExpMed_Masetti"      = "Masetti2022",
  "P21_2024_NPJPrecisOncol_Cheng" = "Cheng2024",
  "P08_2020_CommunBiol_Dong"      = "Dong2020",
  "P09_2022_EurUrol_Cheng"        = "Cheng2022",
  "P12_2021_CancerCell_Kfoury"    = "Kfoury2021",
  "P07_2022_Science_Chan"         = "Chan2022"
)

study_col <- setNames(
  c("#f58231","#bcf60c","#46f0f0","#911eb4","#999999",
    "#f032e6","#4363d8","#2b1700","#e6194B","#ffe119",
    "#009CD2","#3cb44b","#008080","#fabebe","#FFFF33",
    "#9A6324","#e6beff"),
  paperorder
)

paper_col <- c(
  "P01"="#009CD2","P02"="#e6194B","P03"="#3cb44b","P05"="#ffe119","P06"="#4363d8",
  "P07"="#f58231","P08"="#911eb4","P09"="#46f0f0","P10"="#2b1700","P11"="#f032e6",
  "P12"="#bcf60c","P14"="#fabebe","P15"="#008080","P16"="#e6beff","P17"="#9A6324",
  "P19"="#FFFF33","P21"="#999999"
)

etsgroup_col <- c(
  "N"      = "#71A45E", "Benign" = "#87CEFA", "BPH"   = "#0047ab",
  "nonETS" = "#F781BF", "ERG"    = "#984EA3", "ETV"   = "#FF7F00",
  "AR"     = "#E41A1C", "NEPC"   = "#A65628", "DN"    = "#4DAF4A",
  "No_Group" = "#999999"
)

GleasonGrade_col <- c(
  "N"="#71A45E","B"="#87CEFA","BPH"="#0047ab",
  "PCa"="#FFEDA0","GG1"="#FFEDA0","GG2"="#FEB24C","GG2_GG3"="#FEB24C",
  "GG3"="#FC4E2A","GG4"="#E31A1C","GG5"="#800026",
  "mCRPC"="#9E9AC8","mCRPC-NETD"="#6A51A3","mCRPC-teSCC"="#4A1486"
)

Group_col <- c(
  "N"="#71A45E","B"="#87CEFA","BPH"="#0047ab",
  "PCa"="#FD8D3C","mHSPC"="darkgrey","mCRPC"="#6A51A3"
)

group_detailed_col <- c(
  "Luminal"="#E31A1C","Club"="#E2ADB2","Basal"="#1976D2",
  "Hillock"="#87CEFA","Ionocyte"="#4A1486",
  "T cell"="#3B7A57","B cell"="#4DAF4A","Plasma"="#008B8B","NK cell"="#00F0A8",
  "Granulocytes"="#8FBC8F","Macrophage"="#00AB66","Monocyte"="#CEFF00",
  "Neutrophils"="#997950","Dendritic"="#5D4037",
  "Endothelial cell"="#FFD800","SMC"="#FD8D3C","Fibroblast"="#F0E68C",
  "Schwann"="#696969","Adipocytes"="#DEB887"
)

# Extended color map for pseudobulk plots (includes merged group labels)
comp_col <- c(
  etsgroup_col,
  "ETS" = "#984EA3",  # merged ERG+ETV
  "PCa" = "#E41A1C"   # merged all PCa
)

# =============================================================================
# Section 1 — QC Figures (Supplementary Figure 1)
# Source: QC_metrics_raw_prefilter.csv and qc_summary_forsupplementary.csv
# These CSVs are created by 02_data_merge_and_filter.py (Python, Part A) and
# 01_QC_and_annotation.R (R, Part A — scDblFinder doublet counts).
# =============================================================================

message("\n=== Section 1: QC Figures ===")

BASE_DIR <- "/mctp/share/users/hanbyul/projects/human_prostate_sc/cellranger_7.1/17datasets"
PY_FINAL <- file.path(BASE_DIR, "python_final")

# ── 1a. Violin plots: QC metrics per study (pre-filter) ──────────────────────
# Shows the distribution of nCount_RNA, nFeature_RNA, and pct_mt before any
# cell filtering. Red dashed lines = QC thresholds (nCount ≥ 500, nFeature ≥ 250, mt < 20%).

qc_cell <- read.csv(file.path(PY_FINAL, "QC_metrics_raw_prefilter.csv"), row.names = 1)
qc_cell$Paper_ID    <- factor(qc_cell$Paper_ID, levels = paperorder)
qc_cell$study_label <- factor(study_labels[as.character(qc_cell$Paper_ID)],
                               levels = study_labels[paperorder])

make_vln <- function(df, yvar, ylab, fill_col,
                     threshold = NULL, log10_scale = FALSE,
                     keep_direction = "above") {
  p <- ggplot(df, aes(x = study_label, y = .data[[yvar]], fill = Paper_ID)) +
    geom_violin(scale = "width", trim = TRUE, linewidth = 0.3) +
    stat_summary(fun = median, geom = "crossbar", width = 0.5,
                 linewidth = 0.5, color = "black", fatten = 1) +
    scale_fill_manual(values = fill_col) +
    labs(x = NULL, y = ylab) +
    theme_classic(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
          legend.position = "none",
          panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3))

  if (log10_scale) p <- p + scale_y_log10()

  if (!is.null(threshold)) {
    rect_ymin <- if (log10_scale && keep_direction == "above") threshold / 10 else -Inf
    rect_ymax <- if (keep_direction == "above") threshold else Inf
    p <- p +
      annotate("rect", xmin = -Inf, xmax = Inf,
               ymin = rect_ymin, ymax = rect_ymax,
               fill = "red", alpha = 0.08) +
      geom_hline(yintercept = threshold, linetype = "dashed",
                 color = "red", linewidth = 0.6)
    if (log10_scale && keep_direction == "above")
      p <- p + expand_limits(y = threshold / 10)
  }
  p
}

p_count   <- make_vln(qc_cell, "total_counts",      "nCount_RNA (log10)",   study_col,
                       threshold = 500,  log10_scale = TRUE,  keep_direction = "above")
p_feature <- make_vln(qc_cell, "n_genes_by_counts", "nFeature_RNA (log10)", study_col,
                       threshold = 250,  log10_scale = TRUE,  keep_direction = "above")
p_mt      <- make_vln(qc_cell, "pct_counts_mt",     "Percent.mt",           study_col,
                       threshold = 20,   log10_scale = FALSE, keep_direction = "below")

pdf(file.path(PY_FINAL, "SFig1a_QC_violin_prefilter.pdf"), width = 14, height = 10)
gridExtra::grid.arrange(p_count, p_feature, p_mt, ncol = 1)
dev.off()
message("  Saved: SFig1a_QC_violin_prefilter.pdf")

# ── 1b. Bar chart: cell counts at each filtering stage per study ──────────────
# Four stages: Before QC, After QC (thresholds), Singlets (doublets removed),
# After annotation QC (lowqual/unknown cell types removed).

qc_sum <- read.csv(file.path(PY_FINAL, "qc_summary_forsupplementary.csv"), row.names = 1)
qc_sum <- qc_sum[qc_sum$Exclude %in% "Include ", ]
qc_sum[, c("Cells_After_QC", "Cells_Singlet", "Cells_before_QC")] <-
  lapply(qc_sum[, c("Cells_After_QC", "Cells_Singlet", "Cells_before_QC")], as.numeric)
qc_sum$Study <- factor(qc_sum$Study, levels = paperorder)

# Cells remaining after annotation QC (from merged_a17_hit1_scib_50pc2000g.rds)
# table(seu$ID) with only annotated (non-lowqual) cells
cells_after_annot <- c(
  L001=580,  L002=5471, L003=1331, L004=3065, L005=2431, L006=5116,
  L007=2146, L008=4403, L009=2711, L010=5607, L011=3660, L012=8479,
  L013=2066, L014=4937, L015=2462, L016=4017, L017=6204, L018=9408,
  L019=3934, L020=879,  L021=2304, L022=2311, L023=2763, L024=5811,
  L025=5312, L026=7520, L027=4715, L028=6749, L029=6665, L030=5724,
  L031=1281, L032=4461, L033=1780, L034=6425, L035=4006, L036=8361,
  L037=3625, L038=4435, L039=3688, L040=3996, L041=6656, L042=4991,
  L043=7837, L044=1789, L045=7685, L046=5791, L047=8484, L048=6117,
  L049=8665, L050=9057, L051=4305, L052=3995, L053=4645, L054=2336,
  L055=2865, L056=2110, L057=8617, L058=4894, L059=2073, L060=1383,
  L061=4089, L062=4594, L063=4736, L064=3721, L065=4917, L066=8659,
  L067=5666, L068=5517, L069=4123, L070=5097, L072=5408, L073=3951,
  L074=4530, L075=261,  L076=2536, L077=1898, L078=18897,L079=2012,
  L081=3250, L082=8286, L083=7588, L084=3053, L085=2634, L086=4872,
  L087=3030, L088=3486, L089=5965, L090=3618, L091=933,  L092=3386,
  L093=3251, L094=2409, L095=2478, L096=4026, L097=1242, L098=3573,
  L099=2615, L100=3884, L101=2059, L102=2190, L103=10324,L104=5125,
  L105=2188, L106=3986, L107=1645, L108=2156, L109=1850, L110=3672,
  L111=3035, L112=2255, L113=5112, L114=5605, L115=5534, L116=8312,
  L117=7932, L118=7385, L119=7107, L120=14591,L121=11129,L122=7646,
  L123=6615, L124=8176, L125=9864, L127=5866, L128=7121, L129=6390,
  L130=5927, L131=6243, L132=2128, L133=4146, L134=11199,L135=12296,
  L136=11253,L137=9266, L138=3501, L139=2513, L141=6379, L142=5121,
  L143=4995, L144=3967, L145=4862, L146=6928, L147=5238, L148=11282,
  L149=4756, L150=3461, L151=7270, L152=4918, L153=797,  L154=123,
  L155=1219, L156=963,  L157=447,  L158=319,  L161=76,   L162=750,
  L163=456,  L164=906,  L165=91,   L166=373,  L167=2714, L168=2077,
  L172=1040
)
qc_sum$Cells_after_annotation_QC <- cells_after_annot[qc_sum$ID]

study_bar <- qc_sum %>%
  group_by(Study) %>%
  summarise(
    Before_QC      = sum(Cells_before_QC,          na.rm = TRUE),
    After_QC       = sum(Cells_After_QC,            na.rm = TRUE),
    Singlet        = sum(Cells_Singlet,              na.rm = TRUE),
    After_Annot_QC = sum(Cells_after_annotation_QC, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(study_label = factor(study_labels[as.character(Study)],
                              levels = study_labels[paperorder]))

study_bar_long <- study_bar %>%
  pivot_longer(cols = c(Before_QC, After_QC, Singlet, After_Annot_QC),
               names_to = "Stage", values_to = "Cells") %>%
  mutate(Stage = factor(Stage,
                        levels = c("Before_QC","After_QC","Singlet","After_Annot_QC")))

p_bar <- ggplot(study_bar_long, aes(x = study_label, y = Cells, fill = Stage)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  scale_fill_manual(
    values = c("Before_QC"     = "#FC8D59", "After_QC"      = "#FEE08B",
               "Singlet"       = "#D9D9D9", "After_Annot_QC"= "#006400"),
    labels = c("Before QC", "After QC", "Singlets only", "After annotation QC")
  ) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = NULL, y = "Number of cells", fill = NULL) +
  theme_classic(base_size = 11) +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1, size = 9),
        legend.position = "top",
        panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3))

pdf(file.path(PY_FINAL, "SFig1b_QC_cellcount_bar.pdf"), width = 12, height = 5)
print(p_bar)
dev.off()
message("  Saved: SFig1b_QC_cellcount_bar.pdf")

# ── 1c. Stacked bar: cells removed at each step ───────────────────────────────
study_bar2 <- study_bar %>%
  mutate(
    Removed_QC     = Before_QC   - After_QC,
    Doublets       = After_QC    - Singlet,
    Low_qual_cells = Singlet     - After_Annot_QC,
    Kept           = After_Annot_QC
  ) %>%
  select(study_label, Removed_QC, Doublets, Low_qual_cells, Kept) %>%
  pivot_longer(cols = c(Removed_QC, Doublets, Low_qual_cells, Kept),
               names_to = "Stage", values_to = "Cells") %>%
  mutate(Stage = factor(Stage,
                        levels = c("Removed_QC","Doublets","Low_qual_cells","Kept")))

p_stack <- ggplot(study_bar2, aes(x = study_label, y = Cells, fill = Stage)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  scale_fill_manual(
    values = c("Removed_QC"     = "#FC8D59", "Doublets"       = "#FEE08B",
               "Low_qual_cells" = "#D9D9D9", "Kept"           = "#006400"),
    labels = c("Removed (low QC)", "Doublets removed",
               "Low-quality cell types removed", "Kept (after annotation QC)")
  ) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = NULL, y = "Number of cells", fill = NULL) +
  theme_classic(base_size = 11) +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1, size = 11),
        legend.position = "top",
        legend.text     = element_text(size = 13),
        legend.key.size = unit(0.6, "cm"),
        panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3))

pdf(file.path(PY_FINAL, "SFig1c_QC_cellcount_stacked.pdf"), width = 12, height = 5)
print(p_stack)
dev.off()
message("  Saved: SFig1c_QC_cellcount_stacked.pdf")

# ── 1d. Doublet rate per sample ───────────────────────────────────────────────
doublet_sample <- qc_sum %>%
  mutate(Doublet_Rate_percent = as.numeric(Doublet_Rate_percent)) %>%
  filter(!is.na(Doublet_Rate_percent)) %>%
  mutate(study_label = factor(study_labels[as.character(Study)],
                              levels = study_labels[paperorder]))

p_doublet <- ggplot(doublet_sample, aes(x = study_label, y = Doublet_Rate_percent)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, linewidth = 0.4, fill = "white") +
  geom_jitter(aes(color = Study), width = 0.15, height = 0, size = 1.8, alpha = 0.9) +
  scale_color_manual(values = study_col) +
  labs(x = NULL, y = "Doublet rate per sample (%)") +
  theme_classic(base_size = 11) +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1, size = 9),
        legend.position = "none",
        panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3))

pdf(file.path(PY_FINAL, "SFig1d_QC_doublet_rate.pdf"), width = 12, height = 4)
print(p_doublet)
dev.off()
message("  Saved: SFig1d_QC_doublet_rate.pdf")

# =============================================================================
# Section 2 — All-cells scANVI UMAP (Figure 1, Supplementary Figure 1b–c)
# Source: merged_a17_hit1_scib_50pc2000g.rds (from step 06)
# scANVI UMAP stored as reduction "umap_scanvi"
# =============================================================================

message("\n=== Section 2: All-cells UMAPs ===")

seu_all <- readRDS(file.path(PY_FINAL, "merged_a17_hit1_scib_50pc2000g.rds"))
scanvi_red_all <- "umap_scanvi"

message("  Cells: ", ncol(seu_all), "  Reduction: ", scanvi_red_all)

# Cell type UMAP (Figure 1b)
pdf(file.path(PY_FINAL, "Fig1_allcells_celltype_umap.pdf"), width = 9, height = 7)
p <- DimPlot(seu_all, reduction = scanvi_red_all,
             group.by = "final_hit1_group_detailed",
             label = TRUE, label.size = 3, repel = TRUE, raster = FALSE) +
  scale_color_manual(values = group_detailed_col, na.value = "lightgrey") +
  ggtitle("All cells — cell type (scANVI)") +
  theme(legend.text = element_text(size = 8), legend.key.size = unit(0.4, "cm"))
print(p)
dev.off()
message("  Saved: Fig1_allcells_celltype_umap.pdf")

# Study UMAP (Supplementary Figure 1b)
pdf(file.path(PY_FINAL, "SFig1b_allcells_study_umap.pdf"), width = 9, height = 7)
p <- DimPlot(seu_all, reduction = scanvi_red_all,
             group.by = "Paper_ID_short",
             label = FALSE, raster = FALSE) +
  scale_color_manual(values = paper_col, na.value = "lightgrey") +
  ggtitle("All cells — study (scANVI)") +
  theme(legend.text = element_text(size = 8), legend.key.size = unit(0.4, "cm"))
print(p)
dev.off()
message("  Saved: SFig1b_allcells_study_umap.pdf")

# Disease group UMAP (Supplementary Figure 1c)
pdf(file.path(PY_FINAL, "SFig1c_allcells_group_umap.pdf"), width = 7, height = 6)
p <- DimPlot(seu_all, reduction = scanvi_red_all,
             group.by = "Group",
             label = FALSE, raster = FALSE) +
  scale_color_manual(values = Group_col, na.value = "lightgrey") +
  ggtitle("All cells — disease group (scANVI)") +
  theme(legend.text = element_text(size = 9), legend.key.size = unit(0.4, "cm"))
print(p)
dev.off()
message("  Saved: SFig1c_allcells_group_umap.pdf")

rm(seu_all); gc()

# =============================================================================
# Section 3 — Luminal cells scANVI UMAP (Figure 2)
# Source: seu_lum_v1_noNoGroup_scib_30pc3000g.rds (from step 06)
# scANVI UMAP stored as reduction "umap_scanvi"
# ETS subgroups: Benign, BPH, ERG, ETV, nonETS, AR, NEPC, DN
# =============================================================================

message("\n=== Section 3: Luminal UMAPs ===")

SCIB_DIR <- file.path(PY_FINAL, "032025_scib")
seu_lum  <- readRDS(file.path(SCIB_DIR, "seu_lum_v1_noNoGroup_scib_30pc3000g.rds"))
scanvi_red_lum <- "umap_scanvi"

message("  Cells: ", ncol(seu_lum), "  Reduction: ", scanvi_red_lum)

# ETS subgroup UMAP (Figure 2a)
pdf(file.path(SCIB_DIR, "Fig2a_lum_etsgroup_umap.pdf"), width = 8, height = 6)
p <- DimPlot(seu_lum, reduction = scanvi_red_lum,
             group.by = "etsgroup",
             label = TRUE, label.size = 3.5, repel = TRUE, raster = FALSE) +
  scale_color_manual(values = etsgroup_col, na.value = "lightgrey") +
  ggtitle("Luminal cells — ETS subgroup (scANVI)") +
  theme(legend.text = element_text(size = 9), legend.key.size = unit(0.4, "cm"))
print(p)
dev.off()
message("  Saved: Fig2a_lum_etsgroup_umap.pdf")

# GleasonGrade UMAP (Figure 2b)
pdf(file.path(SCIB_DIR, "Fig2b_lum_gleason_umap.pdf"), width = 8, height = 6)
p <- DimPlot(seu_lum, reduction = scanvi_red_lum,
             group.by = "GleasonGrade",
             label = FALSE, raster = FALSE) +
  scale_color_manual(values = GleasonGrade_col, na.value = "lightgrey") +
  ggtitle("Luminal cells — Gleason Grade (scANVI)") +
  theme(legend.text = element_text(size = 9), legend.key.size = unit(0.4, "cm"))
print(p)
dev.off()
message("  Saved: Fig2b_lum_gleason_umap.pdf")

# Study UMAP (Supplementary)
pdf(file.path(SCIB_DIR, "SFig_lum_study_umap.pdf"), width = 8, height = 6)
p <- DimPlot(seu_lum, reduction = scanvi_red_lum,
             group.by = "Paper_ID_short",
             label = FALSE, raster = FALSE) +
  scale_color_manual(values = paper_col, na.value = "lightgrey") +
  ggtitle("Luminal cells — study (scANVI)") +
  theme(legend.text = element_text(size = 8), legend.key.size = unit(0.4, "cm"))
print(p)
dev.off()
message("  Saved: SFig_lum_study_umap.pdf")

# =============================================================================
# Section 4 — Feature plots: marker gene / signature expression on luminal UMAP
#             (Supplementary Figure 4d)
# Source: seu_lum_v1_noNoGroup_scib_30pc3000g.rds (already loaded above)
#
# Plots log-normalized expression of known ETS/AR/NE marker genes and
# AR-TCGA / NE-Beltran module scores on the scANVI UMAP.
# Each plot = one gene, colored by quantile-capped expression (10th–90th).
# =============================================================================

message("\n=== Section 4: Luminal feature plots ===")

# Gene sets
AR_TCGA_genes  <- c("KLK3","KLK2","PMEPA1","ABCC4","NKX3-1","C1orf116",
                     "FKBP5","ACSL3","ZBTB10","HERC3","PTGER4","MPHOSPH9",
                     "EAF2","MED28","NNMT","MAF","GNMT","CENPN","ELL2","TMPRSS2")
NE_Beltran_genes <- c("BRINP1","C7orf76","GNAO1","KCNB2","KCND2","LRRC16B","MAP10",
                       "NRSN1","PCSK1","PROX1","RGS7","SCG3","SEC11C","SEZ6",
                       "ST8SIA3","SVOP","SYT11","AURKA","DNMT1","EZH2","MYCN")
single_genes    <- c("AR","ERG","ETV1","ETV4","SPINK1",
                     "AMACR","PCA3","DLX1","HOXC6","SCHLAP1","FOLH1","EZH2","MKI67")

DefaultAssay(seu_lum) <- "RNA"
seu_lum <- NormalizeData(seu_lum, normalization.method = "LogNormalize",
                         scale.factor = 1e4, verbose = FALSE)

genes_present <- rownames(seu_lum)
AR_use <- AR_TCGA_genes[AR_TCGA_genes   %in% genes_present]
NE_use <- NE_Beltran_genes[NE_Beltran_genes %in% genes_present]

seu_lum <- AddModuleScore(seu_lum, features = list(AR_use), name = "AR_TCGA_score",   nbin = 24, seed = 1)
seu_lum <- AddModuleScore(seu_lum, features = list(NE_use), name = "NE_Beltran_score", nbin = 24, seed = 1)
seu_lum$AR_TCGA_score   <- seu_lum$AR_TCGA_score1;   seu_lum$AR_TCGA_score1   <- NULL
seu_lum$NE_Beltran_score <- seu_lum$NE_Beltran_score1; seu_lum$NE_Beltran_score1 <- NULL

gene_features <- single_genes[single_genes %in% genes_present]
all_features  <- c("AR_TCGA_score", "NE_Beltran_score", gene_features)

# Pad UMAP range for consistent axis limits
pad_range <- function(x, frac = 0.03) {
  rx <- range(x, na.rm = TRUE); d <- diff(rx); if (d == 0) d <- 1
  c(rx[1] - d * frac, rx[2] + d * frac)
}

theme_feat <- theme_classic(base_size = 10) +
  theme(aspect.ratio = 1,
        axis.text  = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 12, color = "black"),
        plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
        legend.text  = element_text(size = 9),
        legend.title = element_text(size = 9))

make_feature_plot <- function(seu, reduction, feat) {
  emb  <- Embeddings(seu, reduction)
  xlim <- pad_range(emb[, 1])
  ylim <- pad_range(emb[, 2])

  if (feat %in% colnames(seu@meta.data)) {
    vals <- seu@meta.data[[feat]]
  } else if (feat %in% rownames(seu)) {
    vals <- as.numeric(GetAssayData(seu, slot = "data")[feat, ])
  } else {
    warning(feat, " not found — skipping"); return(NULL)
  }

  lo <- quantile(vals, 0.10, na.rm = TRUE)
  hi <- quantile(vals, 0.90, na.rm = TRUE)
  if (lo >= hi) { lo <- min(vals, na.rm = TRUE); hi <- max(vals, na.rm = TRUE) }
  if (lo >= hi) { lo <- 0; hi <- 1 }

  df <- data.frame(x = emb[, 1], y = emb[, 2], val = vals)
  set.seed(42); df <- df[sample.int(nrow(df)), , drop = FALSE]

  ggplot(df, aes(x, y, colour = val)) +
    geom_point(size = 0.3, shape = 16, stroke = 0) +
    scale_colour_gradientn(
      colours = c("lightgrey","#FEE08B","#F46D43","#A50026"),
      limits  = c(lo, hi), oob = scales::squish, name = feat
    ) +
    coord_fixed(ratio = 1, xlim = xlim, ylim = ylim, clip = "on") +
    labs(x = colnames(emb)[1], y = colnames(emb)[2]) +
    ggtitle(feat) + theme_feat
}

pdf(file.path(SCIB_DIR, "SFig4d_lum_featureplots.pdf"), width = 5, height = 5)
for (feat in all_features) {
  tryCatch({
    p <- make_feature_plot(seu_lum, scanvi_red_lum, feat)
    if (!is.null(p)) print(p)
  }, error = function(e) message("  ERROR — ", feat, ": ", e$message))
}
dev.off()
message("  Saved: SFig4d_lum_featureplots.pdf")

rm(seu_lum); gc()

# =============================================================================
# Section 5 — Pseudobulk bubble (dot) plots (Figure 4c–d)
#
# Uses pseudobulk-level Seurat object where each "cell" = one pseudobulk
# sample from a given donor × ETS subgroup. FindAllMarkers identifies DE genes
# between subgroups; DotPlot visualizes expression × fraction expressed.
#
# Comparisons:
#   Primary PCa: Benign / nonETS / ERG / ETV (Figure 4c)
#   mCRPC:       AR / NEPC / DN              (Figure 4d)
#
# Source: external pseudobulk RDS (generated separately from the single-cell
# data using Seurat's AggregateExpression or a custom pseudobulk pipeline).
# =============================================================================

message("\n=== Section 5: Pseudobulk bubble plots ===")

PSEUDO_RDS <- paste0(
  "/mctp/share/users/hanbyul/projects/tm/human_prostate_sc/",
  "lum_pseudolevel/correctmetapca10_new/seurat_all_allgroups_pca10_0320.rds"
)

if (!file.exists(PSEUDO_RDS)) {
  warning("Pseudobulk RDS not found: ", PSEUDO_RDS, " — skipping Section 5")
} else {
  p1 <- readRDS(PSEUDO_RDS)
  p1$etsgroup <- p1$group  # harmonize column name

  message("  Loaded pseudobulk: ", ncol(p1), " samples")

  # Known marker gene sets for each ETS subgroup
  Benign          <- c("RLN1","RLN2","VGLL3","NKX3-1","MSMB","AZGP1","GSTP1")
  Tumor           <- c("AMACR","PCA3","DLX1","HOXC6","SCHLAP1","FOLH1","EZH2","MKI67")
  known_ETS_pos   <- c("ERG","TDRD1","DLX1","HOXC6","CRISP3","SOX9","EZH2",
                        "MMP9","TFF3","CXCR4","GRPR","ETV1","ETV4","ETV5")
  known_ETS_neg   <- c("SPINK1","ANPEP","RLN1","TFF3","EGFR",
                        "KLK2","KLK3","FKBP5","NKX3-1","ACSL3")
  spatial_ETV     <- c("NEFH","ETV1","KCNN4","SLC45A2","PLAT","ETV5","TWIST1","GRPR","FGF12","ETV4")
  spatial_SPOP    <- c("SRARP","SHISA6","HSD17B6","GLB1L3","ALOX15B","SGCD","AZGP1","PHGR1",
                        "FBXO27","ADAMTSL1","MT1G","TMEM184A","CRACR2B","TNFRSF19","FAM3B",
                        "ARX","SOX14","CLDN8","SLC43A1","DACH2")
  spatial_FOXA1   <- c("RLN1","SYCE1","NLRP9","PNLDC1","CEACAM20","HPGD","ARHGAP36",
                        "ADGRG7","FBXL16","SLC51A","AATK","MUC13","TNNT1","RLN2",
                        "SPESP1","GLYATL1","CHRNA7","CDX2","SLCO1A2","CA13")

  known_markers <- unique(c(
    Benign, Tumor, known_ETS_pos, known_ETS_neg,
    spatial_ETV, spatial_SPOP, spatial_FOXA1,
    AR_TCGA_genes, NE_Beltran_genes,
    "AR","ERG","TDRD1","CRISP3","ETV1","ETV4","ETV5","GRPR",
    "SYP","CHGA","CHGB","AURKA","MYCN","DNMT1","ENO2",
    "RLN1","RLN2","MSMB","AZGP1","GSTP1","VGLL3"
  ))

  group_known_sets <- list(
    Benign = Benign,
    ERG    = c(known_ETS_pos),
    ETV    = c(known_ETS_pos, spatial_ETV),
    nonETS = c(known_ETS_neg, spatial_SPOP, spatial_FOXA1),
    AR     = AR_TCGA_genes,
    NEPC   = NE_Beltran_genes,
    DN     = NE_Beltran_genes
  )

  # Regex to exclude poorly-annotated / lncRNA gene names
  bad_gene_pattern <- "^(AC|AL|AP|BX|CR|CT|CU|FP|HG|KB|KI|LINC|LOC|MIR|MTCO|MTND|RF|RP|Z)[0-9]|\\.[0-9]+$"

  # Smart top-gene picker: hard DE picks + known biology fills
  pick_top <- function(markers_grp, known_genes, n_top, n_pool = 75, n_total) {
    markers_grp <- markers_grp[!grepl(bad_gene_pattern, markers_grp$gene), ]
    markers_grp <- markers_grp[order(-markers_grp$avg_log2FC), ]
    top_de      <- head(markers_grp$gene, n_top)
    pool        <- head(markers_grp$gene, n_pool)
    from_known  <- setdiff(intersect(known_genes, pool), top_de)
    combined    <- unique(c(top_de, from_known))
    combined    <- head(combined, n_total)
    if (length(combined) < n_total) {
      extra    <- markers_grp$gene[!markers_grp$gene %in% combined]
      combined <- c(combined, head(extra, n_total - length(combined)))
    }
    combined
  }

  # Core function: FindAllMarkers → DotPlot variants → save PDFs
  make_bubble_plots <- function(seu, group_col, prefix,
                                colmap           = comp_col,
                                known_markers_vec = known_markers,
                                grp_known_sets   = group_known_sets) {

    Idents(seu) <- seu@meta.data[[group_col]]
    groups_use  <- levels(seu@meta.data[[group_col]])

    message("  >>> FindAllMarkers for: ", prefix)
    markers <- FindAllMarkers(seu, only.pos = TRUE,
                              min.pct = 0.10, logfc.threshold = 0.25, verbose = FALSE)
    if (nrow(markers) == 0) {
      message("  No markers found for ", prefix); return(invisible(NULL))
    }

    gene_cluster <- markers %>%
      arrange(desc(avg_log2FC)) %>%
      distinct(gene, .keep_all = TRUE) %>%
      select(gene, cluster)

    get_top <- function(n_total) {
      n_top <- max(3, floor(n_total * 0.5))
      per_group <- lapply(groups_use, function(grp) {
        sub <- markers[markers$cluster == grp, , drop = FALSE]
        if (nrow(sub) == 0) return(character(0))
        known <- if (grp %in% names(grp_known_sets)) grp_known_sets[[grp]] else character(0)
        pick_top(sub, known_genes = known, n_top = n_top, n_pool = 75, n_total = n_total)
      })
      seen <- character(0); genes_ordered <- character(0)
      for (gv in per_group) {
        new_g <- gv[!gv %in% seen]
        genes_ordered <- c(genes_ordered, new_g)
        seen <- c(seen, new_g)
      }
      genes_ordered
    }

    filter_pct <- function(genes, threshold) {
      if (threshold == 0 || length(genes) == 0) return(genes)
      dp <- DotPlot(seu, features = genes, group.by = group_col)$data
      keep <- dp %>%
        group_by(features.plot) %>%
        summarise(max_pct = max(pct.exp), .groups = "drop") %>%
        filter(max_pct >= threshold) %>%
        pull(features.plot)
      genes[genes %in% as.character(keep)]
    }

    make_label_colors <- function(genes, use_color) {
      if (!use_color) return(rep("black", length(genes)))
      sapply(genes, function(g) {
        cl <- as.character(gene_cluster$cluster[gene_cluster$gene == g][1])
        if (!is.na(cl) && cl %in% names(colmap)) colmap[[cl]] else "grey50"
      })
    }

    make_label_faces <- function(genes, use_color) {
      if (!use_color) return(rep("plain", length(genes)))
      ifelse(genes %in% known_markers_vec, "bold", "plain")
    }

    variants <- list(
      v1 = list(n_total = 10, color = FALSE, pct = 0,  tag = "top10"),
      v2 = list(n_total = 10, color = TRUE,  pct = 0,  tag = "top10_color"),
      v3 = list(n_total = 15, color = TRUE,  pct = 0,  tag = "top15_color"),
      v4 = list(n_total = 15, color = TRUE,  pct = 25, tag = "top15_color_pct25"),
      v5 = list(n_total = 15, color = TRUE,  pct = 50, tag = "top15_color_pct50"),
      v6 = list(n_total = 15, color = TRUE,  pct = 60, tag = "top15_color_pct60")
    )

    flavors <- list(
      list(do_scale = FALSE, tag = "noscale",
           col_scale = scale_color_gradient(low = "grey92", high = "#2500ff",
                                            name = "Avg Exp")),
      list(do_scale = TRUE,  tag = "zscore_redblue",
           col_scale = scale_color_gradient2(low = "#3B4CC0", mid = "white",
                                             high = "#B40426", midpoint = 0,
                                             name = "Avg Exp\n(z-score)")),
      list(do_scale = TRUE,  tag = "zscore_greyblue",
           col_scale = scale_color_gradient(low = "grey92", high = "#2500ff",
                                            name = "Avg Exp\n(z-score)"))
    )

    for (vname in names(variants)) {
      v     <- variants[[vname]]
      genes <- get_top(v$n_total)
      if (v$pct > 0) genes <- filter_pct(genes, v$pct)
      if (length(genes) == 0) { message("  No genes after pct filter — skipping ", vname); next }

      label_cols  <- make_label_colors(genes, v$color)
      label_faces <- make_label_faces(genes, v$color)
      pdf_w       <- max(7, length(genes) * 0.42 + 3)

      for (fl in flavors) {
        out_pdf <- file.path(SCIB_DIR, paste0(prefix, "_", v$tag, "_", fl$tag, ".pdf"))
        p <- DotPlot(seu, features = genes, group.by = group_col,
                     dot.scale = 6, scale = fl$do_scale) +
          fl$col_scale +
          scale_size(range = c(0, 6), limits = c(0, 100)) +
          ggtitle(paste0(prefix, "  [", v$tag, " | ", fl$tag, "]")) +
          xlab(NULL) + ylab(NULL) +
          theme_classic(base_size = 11) +
          theme(axis.text.x  = element_text(angle = 45, hjust = 1, size = 11,
                                            colour = label_cols, face = label_faces),
                axis.text.y  = element_text(size = 12),
                plot.title   = element_text(face = "bold", hjust = 0.5, size = 12),
                legend.text  = element_text(size = 9),
                legend.title = element_text(size = 10))
        pdf(out_pdf, width = pdf_w, height = 3.5)
        print(p)
        dev.off()
        message("    Saved: ", basename(out_pdf))
      }
    }
    invisible(markers)
  }

  # ── Primary PCa: Benign / nonETS / ERG / ETV (Figure 4c) ─────────────────
  p1_pca <- subset(p1, subset = Group %in% c("B","BPH","N","PCa"))
  p1_pca <- subset(p1_pca, subset = etsgroup %in% c("Benign","BPH","nonETS","ERG","ETV"))
  p1_pca$etsgroup <- as.character(p1_pca$etsgroup)
  p1_pca$etsgroup[p1_pca$etsgroup == "BPH"] <- "Benign"  # merge BPH with Benign
  p1_pca$etsgroup <- factor(p1_pca$etsgroup,
                             levels = names(etsgroup_col)[names(etsgroup_col) %in%
                                                            unique(p1_pca$etsgroup)])

  make_bubble_plots(p1_pca, "etsgroup", "Fig4c_pseudo_PCa_4grp",
                    colmap = comp_col, known_markers_vec = known_markers,
                    grp_known_sets = group_known_sets)

  # ── mCRPC: AR / NEPC / DN (Figure 4d) ────────────────────────────────────
  p1_mcrpc <- subset(p1, subset = Group %in% c("mCRPC"))
  p1_mcrpc <- subset(p1_mcrpc, subset = etsgroup %in% c("AR","NEPC","DN"))
  p1_mcrpc$etsgroup <- factor(p1_mcrpc$etsgroup,
                               levels = names(etsgroup_col)[names(etsgroup_col) %in%
                                                              unique(p1_mcrpc$etsgroup)])

  make_bubble_plots(p1_mcrpc, "etsgroup", "Fig4d_pseudo_mCRPC",
                    colmap = comp_col, known_markers_vec = known_markers,
                    grp_known_sets = group_known_sets)

  rm(p1, p1_pca, p1_mcrpc); gc()
}

cat("\nAll figures complete.\n")
cat("Output directory (all cells): ", PY_FINAL, "\n")
cat("Output directory (luminal):   ", SCIB_DIR, "\n")
