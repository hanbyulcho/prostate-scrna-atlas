# Code release header: the 17-type map, the marker panel and small helpers, loaded by 02_cluster_annotate.R, 03_recall_qc_classify.R and the recall step 4 scripts.
# The map reads Sheet1 of celltype_colorcode_update_17papers.xlsx (in this folder): columns detail, detail_reduce and group_detailed. Methods: Cell-type annotation.
# helpers.R -- label mapping, marker panel, small utilities. Sourced after 00_env.R (not included).

# ---- 17-type mapping table: Sheet1 of celltype_colorcode_update_17papers.xlsx ----
CELLMATCH_XLSX <- file.path(C17, "celltype_colorcode_update_17papers.xlsx")
read_cellmatch <- function() openxlsx::read.xlsx(CELLMATCH_XLSX)

# map17_original (recorded only, not used for labels):
#   final_hit1 -> cellmatch$group_detailed
#   then Adipocytes -> Endothelial cell, Neuroendocrine -> Luminal, Neutrophils -> Monocyte
map17_original <- function(hit1, cm) {
  g <- cm$group_detailed[match(hit1, cm$detail)]
  g[g %in% "Adipocytes"] <- "Endothelial cell"
  g[g %in% "Neuroendocrine"] <- "Luminal"
  g[g %in% "Neutrophils"] <- "Monocyte"
  g
}
# map17_revised (the labels): group_detailed of the hit1 row, with the edits below
map17_revised <- function(hit1, cm) {
  g <- cm$group_detailed[match(hit1, cm$detail)]
  g[hit1 %in% "mast"] <- "Mast cell"                 # the table maps hit1 'mast' to 'remove'; here it is Mast cell
  g[g %in% "Adipocytes"] <- "Endothelial cell"       # merged into Endothelial cell
  g[g %in% "Neuroendocrine"] <- "Luminal"            # merged into Luminal
  g[g %in% "Neutrophils"] <- "Neutrophil"            # kept as its own type (CL:0000775)
  g[g %in% "Granulocytes"] <- "Mast cell"            # renamed (CL:0000097)
  g[is.na(g) & !is.na(hit1)] <- "unmapped (hit1 not in table)"
  g
}
fine_revised <- function(hit1, cm) {                 # detail_reduce, with the same renames
  f <- cm$detail_reduce[match(hit1, cm$detail)]
  f[f %in% "Granulocytes"] <- "Mast cell"
  f[f %in% "Neutrophils"] <- "Neutrophil"
  f[is.na(f) & !is.na(hit1)] <- "unmapped (hit1 not in table)"
  f
}
# carried label of the first atlas version, in the release vocabulary (Granulocytes -> Mast cell; fine label Neutrophils -> Neutrophil)
carried_revised <- function(group_detailed, fine_detailed) {
  g <- as.character(group_detailed)
  g[g %in% "Granulocytes"] <- "Mast cell"
  g[as.character(fine_detailed) %in% "Neutrophils"] <- "Neutrophil"
  g
}
NOT_A_TYPE <- c("remove", "unmapped (hit1 not in table)")

# ---- marker panel. Sources (PANEL_SOURCE): A an earlier Methods marker text, B the marker list of the first atlas version's annotation, C the Mast cell and Neutrophil markers. Of the 27 source-A genes, 22 are in dotplot_gene_sources.csv in the deposit, with sources. Source 'annotation decision': the mast cell, monocyte and neutrophil markers named when Neutrophil was kept as a type ----
PANEL <- list(
  "Luminal"          = c("KLK3", "KLK2", "MSMB", "ACPP", "PLA2G2A", "TFF1", "AZGP1", "MT1G", "TSPAN8", "CLDN3"),
  "Basal"            = c("KRT5", "TP63", "KRT14", "DST", "KRT15", "TAGLN", "CYR61", "RGCC", "TIMP3", "CAV1", "MYC"),
  "Club"             = c("SCGB1A1", "SCGB3A1", "LCN2", "PIGR", "PSCA", "FCGBP", "RARRES1", "SPINK1", "CD74", "WFDC2"),
  "Hillock"          = c("KRT13", "APOBEC3A", "CSTB", "LYPD3", "SEPRINB1", "GDF15", "DUOXA2", "GPX2", "CXCL17"),
  "Ionocyte"         = c("FOXI1"),
  "Neuroendocrine"   = c("CHGA", "CHGB", "FOXA2", "ASCL1"),
  "T cell"           = c("CD3D", "CD3E", "CD3G", "CD8A"),
  "B cell"           = c("CD79A", "MS4A1"),
  "Plasma"           = c("IGHG1", "DERL3", "JCHAIN"),
  "NK cell"          = c("NKG7", "GNLY"),
  "Mast cell"        = c("TPSAB1", "CPA3", "MS4A2", "KIT"),
  "Monocyte"         = c("CD14", "S100A8", "S100A9", "FCN1"),
  "Neutrophil"       = c("FCGR3B", "CXCR2", "CSF3R"),
  "Macrophage"       = c("C1QA", "C1QB", "CD68", "C1QC"),
  "Dendritic"        = c("HLA-DRA", "IL3RA", "CLEC4C", "TLR7", "TLR9"),
  "Endothelial cell" = c("PECAM1", "VWF", "ENG", "PLVAP"),
  "Fibroblast"       = c("DCN", "LUM", "FBLN1", "COL1A2"),
  "SMC"              = c("ACTA2", "MYH11", "RGS5", "MYL9"),
  "Schwann"          = c("S100B", "PLP1", "MPZ"),
  "Compartment"      = c("EPCAM", "KRT8", "VIM", "PTPRC"))
PANEL_SOURCE <- c(KLK3="A,B", KLK2="A,B", MSMB="B", ACPP="B", PLA2G2A="B", TFF1="B", AZGP1="B", MT1G="B", TSPAN8="B", CLDN3="B",
  KRT5="A,B", TP63="A,B", KRT14="B", DST="B", KRT15="B", TAGLN="B", CYR61="B", RGCC="B", TIMP3="B", CAV1="B", MYC="B",
  SCGB1A1="A", SCGB3A1="A,B", LCN2="B", PIGR="B", PSCA="B", FCGBP="B", RARRES1="B", SPINK1="B", CD74="B", WFDC2="B",
  KRT13="A,B", APOBEC3A="B", CSTB="B", LYPD3="B", SEPRINB1="B (as written)", GDF15="B", DUOXA2="B", GPX2="B", CXCL17="B",
  FOXI1="A", CHGA="B", CHGB="B", FOXA2="B", ASCL1="B", CD3D="A,B", CD3E="A,B", CD3G="B", CD8A="B", CD79A="A,B", MS4A1="A,B",
  IGHG1="A", DERL3="B", JCHAIN="B", NKG7="A,B", GNLY="B", TPSAB1="annotation decision,B,C", CPA3="annotation decision,B,C", MS4A2="C", KIT="C",
  CD14="annotation decision", S100A8="A,B", S100A9="A,B", FCN1="B", FCGR3B="annotation decision,C", CXCR2="annotation decision,C", CSF3R="C",
  C1QA="A", C1QB="A", CD68="B", C1QC="B", "HLA-DRA"="A", IL3RA="B", CLEC4C="B", TLR7="B", TLR9="B",
  PECAM1="A,B", VWF="A,B", ENG="B", PLVAP="B", DCN="A,B", LUM="A,B", FBLN1="B", COL1A2="B", ACTA2="A,B", MYH11="A", RGS5="B", MYL9="B",
  S100B="A,B", PLP1="A,B", MPZ="B", EPCAM="B", KRT8="B", VIM="B", PTPRC="B")
# Source D: the panels of a marker audit, read at run time from its panels.py (not included), long table only
read_audit_panel <- function() {
  f <- file.path(REV, "annotation_audit_v1/scripts/panels.py")
  txt <- readLines(f)
  out <- list(); sect <- NA
  for (ln in txt) {
    if (grepl("^ATLAS = \\{", ln)) sect <- "ATLAS"
    if (grepl("^EXTRA = \\{", ln)) sect <- "EXTRA"
    if (grepl("^COMPARTMENT = \\{", ln)) sect <- "COMPARTMENT"
    if (grepl("^MIN_CO|^ATLAS_COMPARTMENT|^OTHER|^def ", ln)) sect <- NA
    if (!is.na(sect) && grepl("^\\s*'[^']+'\\s*:\\s*\\[", ln)) {
      nm <- sub("^\\s*'([^']+)'.*", "\\1", ln)
      g <- regmatches(sub("^[^\\[]*\\[([^]]*)\\].*", "\\1", ln), gregexpr("'[^']+'", sub("^[^\\[]*\\[([^]]*)\\].*", "\\1", ln)))[[1]]
      out[[paste0("D_", sect, ":", nm)]] <- gsub("'", "", g)
    }
  }
  out
}

# per-cluster % expressing (counts > 0) and mean log-normalized expression for a gene set
cluster_marker_stats <- function(counts, data, cl, genes) {
  present <- genes[genes %in% rownames(counts)]
  cls <- sort(unique(as.character(cl)))
  res <- list()
  for (k in c(cls, "ALL_CELLS")) {
    idx <- if (k == "ALL_CELLS") seq_along(cl) else which(as.character(cl) == k)
    if (length(present)) {
      ct <- counts[present, idx, drop = FALSE]; dt <- data[present, idx, drop = FALSE]
      pct <- Matrix::rowSums(ct > 0) / length(idx) * 100
      mn  <- Matrix::rowMeans(dt)
    } else { pct <- numeric(0); mn <- numeric(0) }
    res[[k]] <- data.frame(cluster = k, gene = genes, in_matrix = genes %in% rownames(counts),
                           pct_expressing = ifelse(genes %in% present, round(pct[genes], 1), NA),
                           mean_lognorm = ifelse(genes %in% present, round(mn[genes], 3), NA), stringsAsFactors = FALSE)
  }
  do.call(rbind, res)
}
fmt_type <- function(st, genes) {   # "KRT5 85%|1.92; TP63 40%|0.61"; absent genes marked
  s <- st[match(genes, st$gene), ]
  paste(ifelse(s$in_matrix, sprintf("%s %s%%|%s", s$gene, formatC(s$pct_expressing, format = "f", digits = 0),
                                    formatC(s$mean_lognorm, format = "f", digits = 2)),
               paste0(genes, " absent")), collapse = "; ")
}
comp_string <- function(x) {        # "Luminal 120; T cell 3" ordered by count
  x <- x[!is.na(x)]
  if (!length(x)) return("")
  t <- sort(table(x), decreasing = TRUE)
  paste(sprintf("%s %d", names(t), as.integer(t)), collapse = "; ")
}
