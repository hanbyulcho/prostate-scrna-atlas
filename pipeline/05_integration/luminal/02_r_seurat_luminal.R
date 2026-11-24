# Code release header: Seurat v5 CCA, RPCA and Harmony integrations (IntegrateLayers) with UMAPs, for one seed. Output: r_embeddings.h5 and r_seeds.json.
# Input: 10x HDF5 counts and an obs CSV of the cells (from a preparation script not included). Layers are split by study (Paper_ID_short).
# Settings: 3,000 variable features and 50 PCs. Neighbors, clusters (resolution 0.5) and UMAPs use dims 1:30. Methods: Data integration and batch correction.
# --seed N sets set.seed(N) and the seeds of RunPCA, RunUMAP and FindClusters. Two seeds inside IntegrateLayers (CCA and RPCA) stay at 42 in both runs (seeds 0 and 1).
# 02_r_seurat_luminal.R - Seurat v5 CCA / RPCA / Harmony for the luminal benchmark.
#
# Integration calls and arguments: split RNA by Paper_ID_short; NormalizeData(); FindVariableFeatures(nfeatures=3000);
# ScaleData(); RunPCA(npcs=50); FindNeighbors/FindClusters(res 0.5)/RunUMAP on dims 1:30 for pca, integrated.cca,
# integrated.rpca, harmony; IntegrateLayers(CCAIntegration / RPCAIntegration / HarmonyIntegration, orig.reduction="pca",
# all other arguments default). The call order is fixed, because the global RNG state that HarmonyIntegration sees depends on it.
#
# Input and output:
#   * input = 10x-HDF5 counts + obs CSV written by 01_prepare_inputs.py (not included), read with Read10X_h5/CreateSeuratObject
#     into an Assay5 with the raw counts and all genes.
#   * embeddings are written to an HDF5 file (cells x dims) for the Python step.
#   * --seed N (required) calls set.seed(N) before each IntegrateLayers.
#
# Seeds:
#   * --seed N is REQUIRED. It reaches every stochastic call of this step:
#       set.seed(N) at the start and before each IntegrateLayers (CCA, RPCA, Harmony: these use the global RNG;
#                   IntegrateLayers exposes no seed argument);
#       RunPCA(seed.use = 42 + N), RunUMAP(seed.use = 42 + N), FindClusters(random.seed = 0 + N).
#     42 / 42 / 0 are the Seurat 5.1.0 defaults, so N = 0 gives the default calls
#     (plus the explicit set.seed(0) calls). Seeds used are written to r_seeds.json.
#   * the diagnostic PDF colors only by the COLOR_KEYS columns present in the obs CSV (the input count matrix may not carry
#     Group/GleasonGrade; a missing column would stop the step after the integrations).
#
# Fixed seeds and safe writes:
#   * --seed does NOT reach two fixed seeds inside IntegrateLayers (Seurat 5.1.0 source, deparsed at run time and
#     written to r_seeds.json 'seurat_internal_seed_source_lines'): CCA -> FindIntegrationAnchors -> RunCCA() without
#     seed.use, i.e. its default seed.use = 42 (set.seed(42)); RPCA -> RunPCA() per layer without seed.use (default 42).
#     Harmony (harmony 1.2.0 RunHarmony.default) has no seed argument and uses the global RNG set by set.seed(N).
#   * outputs written to temp files, synced and renamed; r_embeddings.h5 (the run wrapper, not included, skips this step when it exists) last.
#
# Usage:
#   $RSCRIPT --no-environ 02_r_seurat_luminal.R <prefix_10x.h5> <prefix_obs.csv> <outdir> --seed N [--no-rds] [--no-pdf]
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("usage: 02_r_seurat_luminal.R <10x.h5> <obs.csv> <outdir> --seed N [--no-rds] [--no-pdf]")
in10x <- args[1]; obs_csv <- args[2]; outdir <- args[3]
seed <- NA
if ("--seed" %in% args) seed <- as.integer(args[which(args == "--seed") + 1])
if (is.na(seed)) stop("--seed N is required")
SEEDS <- list(seed = seed, set_seed_start = seed, set_seed_before_IntegrateLayers = seed,
              RunPCA_seed.use = 42L + seed, RunUMAP_seed.use = 42L + seed, FindClusters_random.seed = 0L + seed)
save_rds <- !("--no-rds" %in% args)
make_pdf <- !("--no-pdf" %in% args)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

libs <- strsplit(Sys.getenv("R_LIBS"), ":")[[1]]
.libPaths(c(libs[nzchar(libs)], .libPaths()))
t0 <- Sys.time()
logm <- function(...) message(format(Sys.time(), "%H:%M:%S"), sprintf(" +%.1fm ", as.numeric(difftime(Sys.time(), t0, units = "mins"))), ...)

suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject); library(hdf5r); library(ggplot2); library(patchwork)
})
options(future.globals.maxSize = 4e9)   # future globals up to 4 GB
logm("R ", R.version.string, "; Seurat ", as.character(packageVersion("Seurat")), "; SeuratObject ",
     as.character(packageVersion("SeuratObject")), "; harmony ", as.character(packageVersion("harmony")),
     "; Matrix ", as.character(packageVersion("Matrix")))
logm("libPaths: ", paste(.libPaths(), collapse = " | "))
logm("seed: ", seed, "; ", paste(names(SEEDS), unlist(SEEDS), sep = "=", collapse = " "))
# fixed seeds inside the integration functions, read from the installed source at run time
seed_lines <- function(f) {
  txt <- deparse(eval(parse(text = f)))
  hit <- trimws(grep("set\\.seed|seed\\.use|random\\.seed|future\\.seed|RunPCA\\(|RunCCA\\(", txt, value = TRUE))
  if (length(hit)) hit else "no seed argument, set.seed call or RunPCA/RunCCA call in the deparsed source"
}
SEEDS$seurat_internal_seed_source_lines <- sapply(c("Seurat::CCAIntegration", "Seurat::RPCAIntegration",
  "Seurat::HarmonyIntegration", "Seurat::FindIntegrationAnchors", "Seurat:::RunCCA.default", "Seurat:::RunPCA.default",
  "harmony:::RunHarmony.default"), seed_lines, simplify = FALSE)
SEEDS$seurat_internal_fixed_seeds <- paste(
  "Not changed by --seed: CCA -> FindIntegrationAnchors calls RunCCA() without seed.use (RunCCA.default seed.use = 42,",
  "set.seed(42)); RPCA -> RunPCA() per layer without seed.use (RunPCA.default seed.use = 42). Harmony: no seed",
  "argument; global RNG (set.seed(N) before IntegrateLayers). Seurat", as.character(packageVersion("Seurat")),
  ", harmony", as.character(packageVersion("harmony")))
# thread settings of this step (set by the run scripts, not included), recorded with the seeds
SEEDS$threads_env <- as.list(Sys.getenv(c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "RCPP_PARALLEL_NUM_THREADS",
                                          "THREADS_PER_SEED_RULE", "DRYRUN_THREADS")))
SEEDS$visible_cores <- length(parallel::mcaffinity())
writeLines(jsonlite::toJSON(SEEDS, auto_unbox = TRUE, pretty = TRUE), file.path(outdir, "r_seeds.json"))

BATCH_KEY <- "Paper_ID_short"
source("<BASE_DIR>/17datasets/colorcode.R")  # color palettes (colorcode.R, not included)
COLOR_KEYS <- c("etsgroup", "Group", "GleasonGrade", "Paper_ID_short")
palette_for_key <- list(etsgroup = etsgroup_col, Group = Group_col, GleasonGrade = GleasonGrade_col,
                        Paper_ID_short = paper_col)

# ---- load -------------------------------------------------------------------------------------------------------
counts <- Read10X_h5(in10x, use.names = TRUE, unique.features = TRUE)
seu <- CreateSeuratObject(counts = counts, assay = "RNA")
rm(counts); invisible(gc())
cell_df <- read.csv(obs_csv, row.names = 1, check.names = FALSE)
stopifnot(all(colnames(seu) %in% rownames(cell_df)))
cell_df <- cell_df[colnames(seu), , drop = FALSE]
seu <- AddMetaData(seu, metadata = cell_df)
if (!inherits(seu[["RNA"]], "Assay5")) seu[["RNA"]] <- as(seu[["RNA"]], "Assay5")
logm("cells ", ncol(seu), " genes ", nrow(seu), " batches ", length(unique(seu[[BATCH_KEY, drop = TRUE]])))

# ---- run_integration() body, same calls (seed arguments added; equal to the defaults for seed 0) ---------------
seed_if <- function(tag) { if (!is.na(seed)) set.seed(seed); logm(tag, " RNG digest ", if (exists(".Random.seed")) sum(.Random.seed[1:10] %% 9973) else "none") }
seed_if("start")
seu[["RNA"]] <- split(seu[["RNA"]], f = seu[[BATCH_KEY, drop = TRUE]])
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, nfeatures = 3000)
seu <- ScaleData(seu)
seu <- RunPCA(seu, npcs = 50, seed.use = SEEDS$RunPCA_seed.use)
logm("PCA done")
seu <- FindNeighbors(seu, dims = 1:30, reduction = "pca")
seu <- FindClusters(seu, resolution = 0.5, cluster.name = "unintegrated_clusters", random.seed = SEEDS$FindClusters_random.seed)
seu <- RunUMAP(seu, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated", seed.use = SEEDS$RunUMAP_seed.use)
logm("unintegrated done")

seed_if("CCA")
seu <- IntegrateLayers(seu, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca", verbose = FALSE)
seu <- FindNeighbors(seu, reduction = "integrated.cca", dims = 1:30)
seu <- FindClusters(seu, resolution = 0.5, cluster.name = "cca_clusters", random.seed = SEEDS$FindClusters_random.seed)
seu <- RunUMAP(seu, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca", seed.use = SEEDS$RunUMAP_seed.use)
logm("CCA done")

seed_if("RPCA")
seu <- IntegrateLayers(seu, method = RPCAIntegration, orig.reduction = "pca", new.reduction = "integrated.rpca", verbose = FALSE)
seu <- FindNeighbors(seu, reduction = "integrated.rpca", dims = 1:30)
seu <- FindClusters(seu, resolution = 0.5, cluster.name = "rpca_clusters", random.seed = SEEDS$FindClusters_random.seed)
seu <- RunUMAP(seu, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca", seed.use = SEEDS$RunUMAP_seed.use)
logm("RPCA done")

seed_if("Harmony")
seu <- IntegrateLayers(seu, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony", verbose = FALSE)
seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:30)
seu <- FindClusters(seu, resolution = 0.5, cluster.name = "harmony_clusters", random.seed = SEEDS$FindClusters_random.seed)
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", seed.use = SEEDS$RunUMAP_seed.use)
logm("Harmony done")

COLOR_KEYS <- intersect(COLOR_KEYS, colnames(seu@meta.data))   # only columns present in the obs CSV
if (make_pdf) {
  methods_info <- list(list(reduction = "umap.unintegrated", label = "Unintegrated"), list(reduction = "umap.cca", label = "CCA"),
                       list(reduction = "umap.rpca", label = "RPCA"), list(reduction = "umap.harmony", label = "Harmony"))
  pdf(file.path(outdir, "umap_batch_luminal_R.pdf"), width = 5 * length(COLOR_KEYS), height = 5)
  for (m in methods_info) {
    plots <- lapply(COLOR_KEYS, function(key) {
      cols <- palette_for_key[[key]]
      vals <- as.character(unique(seu[[key, drop = TRUE]]))
      cols <- cols[names(cols) %in% vals]
      p <- DimPlot(seu, reduction = m$reduction, group.by = key, label = FALSE, raster = FALSE) +
        ggtitle(paste0(m$label, " — ", key)) +
        theme(legend.text = element_text(size = 7), legend.key.size = unit(0.4, "cm"))
      if (length(cols) > 0) p <- p + scale_color_manual(values = cols, na.value = "lightgrey")
      p
    })
    print(wrap_plots(plots, nrow = 1))
  }
  dev.off()
}

seu <- JoinLayers(seu)

# ---- outputs (temp file -> sync -> size check -> rename; r_embeddings.h5 last) -----------------------------
finalize <- function(tmp, final) {
  system2("sync", tmp)
  if (!file.exists(tmp) || file.size(tmp) == 0) stop("empty output ", tmp)
  if (!file.rename(tmp, final)) stop("rename failed ", tmp)
}
TMP <- paste0(".tmp", Sys.getpid())
red_map <- c(pca = "X_pca_unintegrated", umap.unintegrated = "X_umap_unintegrated", integrated.cca = "X_cca",
             umap.cca = "X_umap_cca", integrated.rpca = "X_rpca", umap.rpca = "X_umap_rpca", harmony = "X_harmony",
             umap.harmony = "X_umap_harmony")
vf <- VariableFeatures(seu)
writeLines(vf, file.path(outdir, paste0("r_variable_features.txt", TMP)))
finalize(file.path(outdir, paste0("r_variable_features.txt", TMP)), file.path(outdir, "r_variable_features.txt"))
write.csv(seu@meta.data[, c("unintegrated_clusters", "cca_clusters", "rpca_clusters", "harmony_clusters")],
          file.path(outdir, paste0("r_clusters.csv", TMP)))
finalize(file.path(outdir, paste0("r_clusters.csv", TMP)), file.path(outdir, "r_clusters.csv"))
if (save_rds) {
  saveRDS(seu, file.path(outdir, paste0("luminal_seurat_integrated.rds", TMP)))
  finalize(file.path(outdir, paste0("luminal_seurat_integrated.rds", TMP)), file.path(outdir, "luminal_seurat_integrated.rds"))
}
writeLines(capture.output(sessionInfo()), file.path(outdir, "r_sessionInfo.txt"))
h5f <- file.path(outdir, "r_embeddings.h5")
h5t <- paste0(h5f, TMP)
if (file.exists(h5t)) file.remove(h5t)
h <- H5File$new(h5t, mode = "w")
h[["cells"]] <- colnames(seu)
for (r in names(red_map)) {
  e <- Embeddings(seu, reduction = r)
  stopifnot(identical(rownames(e), colnames(seu)))
  h[[red_map[[r]]]] <- t(e)      # hdf5r writes R column-major; transposed so h5py reads cells x dims
}
h$close_all()
finalize(h5t, h5f)
logm("wrote ", h5f, "; variable features ", length(vf), "; rds ", save_rds)
