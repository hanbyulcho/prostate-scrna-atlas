# Code release header: infers arm-level copy-number values of each luminal unit from pseudobulk expression. Methods: Copy-number inference.
# Input: the outputs of 01_pseudobulk_v4.R, sample_metadata.csv, gene positions, hg38 cytobands and the Chan et al. 2022 alteration table.
# Output: arm_values_full_v4.csv.gz, plus locus values, half-split arm values and other tables.
# Per unit: log2(CPM + 1) after TMM, genes above 2 in unit and reference, log ratio clipped to [-3, 3] and median-centered, then DNAcopy (set.seed(1)).
# Arm value: mean segment value over the arm's genes, missing below 20 genes. ref_ids() defines four reference options. The labeling rule uses option b: benign luminal pseudobulks with 50 or more luminal cells and the same chemistry, own donor included.
# 02_cnv_v4.R
# Variant P: TMM log2(CPM+1); genes > 2 in unit and reference; log ratio clipped +-3; median-centered; DNAcopy smooth.CNA + segment (defaults; set.seed(1));
# arm value = gene-weighted mean of segment means (>= 20 genes). All D4 and D3 units under options a-d; all D4 half-splits.
.libPaths(c("<PROJECT_DIR>/redo_env_v1/Rlib_R433_overlay",
            "<R_LIBRARY>", .Library))
suppressPackageStartupMessages({ library(data.table); library(Matrix); library(edgeR); library(DNAcopy); library(parallel) })
args <- commandArgs(TRUE); NCORES <- if (length(args)) as.integer(args[1]) else 40L
t_start <- Sys.time()
W   <- "<input folder: the work folder of 01_pseudobulk_v4.R (pb_counts_v4.rds, pb_meta_v4.csv, halves/)>"
OUT <- "<output folder: arm values, locus values and half-split arm values>"; dir.create(OUT, showWarnings = FALSE)
R   <- "<PROJECT_DIR>"
CPIL <- "<input folder: holds data/chan2022_molecular_alterations.csv, the S19 alteration table of Chan et al. 2022>"
MIN_ARM_GENES <- 20L; LOCUS_WIN <- 3e6; CLIP <- 3; EXPR_MIN <- 2
YOUNG <- c("L126","L127","L128","L129","L130","L131","L132"); HIRZ <- c("L037","L041","L042","L043"); TS <- c("L138","L139")

X <- readRDS(file.path(W, "pb_counts_v4.rds")); META <- fread(file.path(W, "pb_meta_v4.csv")); stopifnot(identical(colnames(X), META$pb_id))
MD <- fread(file.path(R, "repository_files_v28_2026-11-05/sample_metadata.csv"))[, .(library = internal_sample_id, study_short_id, patient_id, disease_state, library_chemistry, pretreatment_adt)]
stopifnot(!anyDuplicated(MD$library))
MD[, chem_group := fifelse(library_chemistry %in% c("10x 3' v3", "10x 3' v3.1"), "10x 3' v3/v3.1", library_chemistry)]
META <- merge(META, MD, by = "library", all.x = TRUE, sort = FALSE); META <- META[match(colnames(X), pb_id)]; stopifnot(identical(colnames(X), META$pb_id), !anyNA(META$disease_state))
META[, pkey := paste(study_short_id, patient_id, sep = ":")]
Gall <- fread(file.path(R, "cnv_pseudobulk_pilot_v1/work/gene_positions_gtf.tsv"), header = FALSE, col.names = c("gene_id","gene_name","chrom","start","end","gene_type"))
G <- Gall[chrom %in% paste0("chr", c(1:22, "X"))][!grepl("^HLA-|^IG[HKL][VDJC]|^RP[LS][0-9]|^MT-", gene_name)]
G[, chr_i := match(chrom, paste0("chr", c(1:22, "X")))]
cb <- as.data.table(readRDS("<R_LIBRARY>/circlize/extdata/cytoband_list.rds")[["hg38"]]); setnames(cb, c("chrom","bstart","bend","band","stain"))
cen <- cb[substr(band, 1, 1) == "q", .(cen = min(bstart)), by = chrom]
G <- merge(G, cen, by = "chrom"); G[, arm := paste0(sub("chr", "", chrom), fifelse(start < cen, "p", "q"))]
setorder(G, chr_i, start); G <- G[gene_id %in% rownames(X)]

## TMM over all full pseudobulks; halves use their parent unit's factor (project choice)
nf <- calcNormFactors(DGEList(counts = X))$samples$norm.factors; names(nf) <- colnames(X)
LC <- sweep(X[G$gene_id, ], 2, colSums(X) * nf, "/") * 1e6; LC <- log2(LC + 1)
cat(format(Sys.time(), "%T"), "TMM done, pseudobulks", ncol(X), "\n")

U <- META[kind %in% c("unit_d4","unit_d3")]
U[, definition := fifelse(kind == "unit_d4", "D4", "D3")]
U[, class := fifelse(library %in% YOUNG, "young donor", fifelse(library %in% HIRZ, "older normal (Hirz)", fifelse(library %in% TS, "older donor (TS)",
             fifelse(disease_state == "benign", "benign", fifelse(disease_state == "BPH", "BPH", fifelse(disease_state == "mCRPC", "mCRPC", fifelse(disease_state == "primary_tumor", "primary tumour", "other normal")))))))]
print(U[, .N, by = .(definition, class)][order(definition, class)]); stopifnot(!any(U$class == "other normal"))
libl <- META[kind == "lib_luminal" & n_cells >= 50]
ref_ids <- function(u, opt) {
  if (opt == "a") return(META[kind == "nonepi_celltype" & library == u$library]$pb_id)
  if (opt == "b") { if (is.na(u$chem_group) || u$chem_group == "") return(character(0)); return(libl[disease_state == "benign" & chem_group == u$chem_group]$pb_id) }
  if (opt == "c") return(libl[library %in% YOUNG & pkey != u$pkey]$pb_id)
  if (opt == "d") return(libl[pkey == u$pkey & library != u$library & disease_state %in% c("benign","normal","BPH")]$pb_id) }
JOBS <- rbindlist(lapply(seq_len(nrow(U)), function(i) rbindlist(lapply(c("a","b","c","d"), function(o) {
  r <- ref_ids(U[i], o); rm_ <- META[match(r, pb_id)]
  data.table(pb_id = U$pb_id[i], option = o, n_ref_groups = length(r), ref_ids = paste(r, collapse = ";"), ref_libraries = paste(unique(rm_$library), collapse = ";"), ref_cells = sum(rm_$n_cells),
             ref_tumour_type_cells = if (o != "a") sum(rm_$n_tumour_type) else NA_real_,
             ref_includes_other_chem_group = length(r) > 0 && any(rm_$chem_group != U$chem_group[i])) }))))
JOBS[, assessable := n_ref_groups > 0]
REF <- lapply(split(JOBS[assessable == TRUE], by = c("pb_id","option")), function(j) { r <- strsplit(j$ref_ids, ";")[[1]]; if (length(r) == 1) LC[, r] else rowMeans(LC[, r, drop = FALSE]) })

## loci: driver loci, 21q22 interval, Chan S19 copy-number rows (an exported copy, checked equal to the xlsx)
drv <- data.table(locus = c("NKX3-1 (8p21.2)","PTEN (10q23.31)","RB1 (13q14.2)","MAP3K7 (6q15)","CHD1 (5q21.1)","MYC (8q24.21)"), gene = c("NKX3-1","PTEN","RB1","MAP3K7","CHD1","MYC"))
drv <- merge(drv, Gall[, .(gene = gene_name, chrom, gstart = start, gend = end)], by = "gene", sort = FALSE)
drv[, `:=`(wstart = gstart - LOCUS_WIN, wend = gend + LOCUS_WIN, exclude = gene)]
erg <- Gall[gene_name == "ERG"]; tmp <- Gall[gene_name == "TMPRSS2"]; I21 <- c(min(erg$start, tmp$start), max(erg$end, tmp$end))
S19 <- fread(file.path(CPIL, "data/chan2022_molecular_alterations.csv")); setnames(S19, c("sample","gene","type","alteration","cdna","value"))
S19 <- S19[type %in% c("Amplification","Gain","Deletion","Loss") & (is.na(alteration) | alteration %in% c("-","") | grepl("^[0-9X]+[pq][0-9.]+$", alteration))]
S19[, direction := fifelse(type %in% c("Amplification","Gain"), "gain", "loss")]
alias <- c(CDKN2Ap16INK4A = "CDKN2A", CDKN2Ap14ARF = "CDKN2A", B2m = "B2M")     # project choice: symbol variants
S19[, gtf_gene := fifelse(gene %in% names(alias), alias[gene], gene)]
gpos <- Gall[!duplicated(gene_name), .(gtf_gene = gene_name, chrom, gstart = start, gend = end)]
S19 <- merge(S19, gpos, by = "gtf_gene", all.x = TRUE, sort = FALSE)
S19[, locus_from := fifelse(!is.na(chrom), "gene +-3 Mb", "reported cytoband")]     # only MAPK31 lacks a genes.gtf symbol
for (i in which(is.na(S19$chrom))) { a <- S19$alteration[i]; ch <- paste0("chr", sub("([0-9X]+)[pq].*", "\\1", a)); bd <- sub("^[0-9X]+", "", a)
  z <- cb[chrom == ch & band == bd]; if (nrow(z)) set(S19, i, c("chrom","gstart","gend"), list(ch, z$bstart[1] + LOCUS_WIN, z$bend[1] - LOCUS_WIN)) }
S19[is.na(chrom), locus_from := "not locatable"]; S19[, cytoband := alteration]
for (i in which(is.na(S19$alteration) | S19$alteration %in% c("-",""))) { z <- cb[chrom == S19$chrom[i] & bstart <= S19$gstart[i] & bend > S19$gstart[i]]
  set(S19, i, "cytoband", paste0(sub("chr", "", S19$chrom[i]), z$band[1], " (not in S19; from genes.gtf position)")) }
cat("S19 copy-number rows:", nrow(S19), " not locatable:", sum(S19$locus_from == "not locatable"), "\n")
LOCI <- unique(S19[locus_from != "not locatable", .(gene, gtf_gene, chrom, wstart = gstart - LOCUS_WIN, wend = gend + LOCUS_WIN, exclude = fifelse(locus_from == "gene +-3 Mb", gtf_gene, ""))])
LOCI <- rbind(drv[, .(gene, gtf_gene = gene, chrom, wstart, wend, exclude)], LOCI)[!duplicated(gene)]
CHANLIB <- unique(META[study_short_id == "P07"]$library)

profile <- function(cnt, eff, vr, with_loci) {
  vu <- log2(cnt[G$gene_id] / eff * 1e6 + 1); keep <- vu > EXPR_MIN & vr > EXPR_MIN; g <- G[keep]
  d <- pmax(pmin(vu[keep] - vr[keep], CLIP), -CLIP); x <- d - median(d)
  cna <- smooth.CNA(CNA(cbind(x), g$chr_i, g$start, data.type = "logratio", sampleid = "u"))
  set.seed(1); sg <- segment(cna, verbose = 0)$output; gseg <- rep(sg$seg.mean, sg$num.mark); stopifnot(length(gseg) == nrow(g))
  arm <- data.table(arm = g$arm, v = gseg)[, .(value = mean(v), n_genes = .N), by = arm][n_genes < MIN_ARM_GENES, value := NA_real_]
  out <- list(arm = arm, n_genes = nrow(g), n_seg = nrow(sg))
  if (with_loci) {
    out$loc <- rbindlist(lapply(seq_len(nrow(LOCI)), function(k) { z <- LOCI[k]
      ii <- which(g$chrom == z$chrom & g$start <= z$wend & g$end >= z$wstart & g$gene_name != z$exclude)
      data.table(gene = z$gene, value = if (length(ii)) mean(gseg[ii]) else NA_real_, n_genes = length(ii)) }))
    in21 <- g$chrom == "chr21" & g$start <= I21[2] & g$end >= I21[1] & !g$gene_name %in% c("ERG","TMPRSS2"); rest <- g$arm == "21q" & !in21 & !g$gene_name %in% c("ERG","TMPRSS2")
    out$c21 <- if (sum(in21) && sum(rest)) mean(x[in21]) - mean(x[rest]) else NA_real_; out$c21_n <- sum(in21) }
  out }

FJ <- JOBS[assessable == TRUE]; t0 <- Sys.time()
full <- mclapply(seq_len(nrow(FJ)), function(j) { jb <- FJ[j]; u <- jb$pb_id; p <- profile(X[, u], sum(X[, u]) * nf[u], REF[[paste(u, jb$option, sep = ".")]], TRUE)
  list(arm = data.table(pb_id = u, option = jb$option, p$arm), loc = data.table(pb_id = u, option = jb$option, p$loc),
       uni = data.table(pb_id = u, option = jb$option, n_genes = p$n_genes, n_segments = p$n_seg, c21_raw_contrast = p$c21, c21_n_genes = p$c21_n)) }, mc.cores = NCORES, mc.preschedule = TRUE)
bad <- which(vapply(full, inherits, logical(1), "try-error")); if (length(bad)) stop("full job failed: ", full[[bad[1]]])
t_full <- as.numeric(difftime(Sys.time(), t0, units = "secs")); cat(format(Sys.time(), "%T"), "full units done", nrow(FJ), "jobs,", round(t_full), "s\n")
ARM <- rbindlist(lapply(full, `[[`, "arm")); LOC <- rbindlist(lapply(full, `[[`, "loc")); UNI <- rbindlist(lapply(full, `[[`, "uni"))
LOC <- LOC[sub("\\|.*", "", pb_id) %in% CHANLIB | gene %in% drv$gene]; rm(full); invisible(gc())

wr <- function(d, name) { f <- file.path(OUT, name); gz <- grepl("\\.gz$", name); tmp <- paste0(f, ".tmp")
  fwrite(d, tmp, compress = if (gz) "gzip" else "none"); system2("sync", tmp); stopifnot(file.size(tmp) > 0, nrow(fread(tmp)) == nrow(d))
  if (gz) stopifnot(identical(readBin(tmp, "raw", 2), as.raw(c(0x1f, 0x8b)))); file.rename(tmp, f) }
info <- U[, .(pb_id, definition, library, source, cluster, n_cells, total_umi, n_tumour_type, n_no_fs_row, class, study_short_id, patient_id, disease_state, library_chemistry, chem_group, pretreatment_adt)]
wr(info, "units_v4.csv"); wr(JOBS, "unit_reference_sets_v4.csv.gz"); wr(ARM, "arm_values_full_v4.csv.gz"); wr(LOC, "locus_values_v4.csv.gz"); wr(UNI, "unit_summary_v4.csv.gz")
wr(S19[, .(sample, gene, gtf_gene, type, alteration, cytoband, value, direction, chrom, locus_from)], "chan_s19_truth_rows_v4.csv")
cat(format(Sys.time(), "%T"), "full tables written\n")

## half-splits, D4 units only: one worker per unit
HIx <- fread(file.path(W, "halves_index_v4.csv"))
FJ4 <- FJ[pb_id %in% U[definition == "D4"]$pb_id]; uj <- unique(FJ4$pb_id); uj <- uj[order(-U$n_cells[match(uj, U$pb_id)])]
t0 <- Sys.time()
halves <- mclapply(uj, function(u) {
  L <- sub("\\|.*", "", u); H <- readRDS(file.path(W, "halves", paste0(L, ".rds"))); cols <- grep(paste0("^", gsub("|", "\\|", u, fixed = TRUE), "\\|s"), colnames(H), value = TRUE)
  H <- H[, cols, drop = FALSE]; opts <- FJ4[pb_id == u]$option; res <- list()
  for (cc in cols) { cnt <- as.numeric(H[, cc]); names(cnt) <- rownames(H); eff <- sum(cnt) * nf[u]
    for (o in opts) { p <- profile(cnt, eff, REF[[paste(u, o, sep = ".")]], FALSE)
      res[[length(res) + 1]] <- data.table(pb_id = u, split = as.integer(sub(".*\\|s([0-9]+)\\|h.*", "\\1", cc)), half = as.integer(sub(".*\\|h", "", cc)), option = o, p$arm) } }
  rbindlist(res) }, mc.cores = NCORES, mc.preschedule = FALSE)
bad <- which(vapply(halves, inherits, logical(1), "try-error")); if (length(bad)) stop("half job failed: ", halves[[bad[1]]])
HARM <- rbindlist(halves); t_half <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(format(Sys.time(), "%T"), "halves done:", nrow(HARM), "arm rows,", round(t_half), "s\n")
wr(HARM, "arm_values_halves_v4.csv.gz")
wr(data.table(step = c("full units (D4 + D3, all options)", "half-splits (D4, all options)"), jobs = c(nrow(FJ), nrow(HARM[, .N, by = .(pb_id, split, half, option)])), wall_s = c(t_full, t_half), cores = NCORES), "runtime_cnv_v4.csv")
cat("total_s", round(as.numeric(difftime(Sys.time(), t_start, units = "secs"))), "\n")
