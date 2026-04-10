"""
05_integration_luminal_cells.py
================================
Step 5 of the pipeline: batch correction and benchmarking for luminal cells.

Cho et al. — A harmonized single-cell RNA-seq atlas of the human prostate
Scientific Data (in submission)

Workflow position:
  03_batch_correction_luminal.R → [THIS FILE] → 06_h5ad_to_seurat.R

Input:
  seu_lum_v1_noNoGroup.h5ad        — exported by 03_batch_correction_luminal.R
    Contains Seurat v5 embeddings from CCA, RPCA, and Harmony integration,
    plus unintegrated PCA. All stored in adata.obsm with Seurat naming:
      X_pca              → renamed to X_pca_unintegrated   (Step 2 below)
      X_umap.unintegrated → renamed to X_umap_unintegrated
      X_integrated.cca   → renamed to X_cca
      X_umap.cca         → renamed to X_umap_cca
      X_integrated.rpca  → renamed to X_rpca
      X_umap.rpca        → renamed to X_umap_rpca
      X_harmony          → kept as X_harmony
      X_umap.harmony     → renamed to X_umap_harmony
  seu_lum_v1_noNoGroup_metadata.csv — full obs metadata (re-attached after R export)

Parameters (manuscript: luminal cells):
  - 3,000 highly variable genes (batch-aware, Seurat v3 flavor, pseudo-counts)
  - 30 principal components
  - Batch key : Paper_ID_short (study short code)
  - Label key : etsgroup (ETS molecular subgroup — used for bio conservation)
  - kNN neighbors for UMAP: 30

Integration methods (8 total):
  From R (03_batch_correction_luminal.R):
    1. Unintegrated (PCA, Seurat)
    2. CCA   — Canonical Correlation Analysis (Seurat IntegrateLayers)
    3. RPCA  — Reciprocal PCA (Seurat IntegrateLayers)
    4. Harmony (Seurat IntegrateLayers)
  From Python (this script):
    5. Scanorama           — panoramic stitching (dimred=30)
    6. BBKNN               — graph-only; UMAP as proxy †
    7. scVI                — n_latent=30, n_layers=2, max_epochs=400, early stopping
    8. scANVI              — n_latent=30, max_epochs=20 ← selected method

† See note in Step 5.

Benchmarking: scib-metrics (Luecken et al. 2022, Nat Methods)
  Label key for bio conservation: etsgroup (ETS molecular subgroup)
  Total = 0.4 × Batch correction + 0.6 × Bio conservation
  Results → scib_results_raw.csv / scib_results_minmax.csv (Supplementary Tables 5, 6)

Output:
  seu_lum_v1_noNoGroup_scib_30pc3000g.h5ad      — all embeddings, X = lognorm
  seu_lum_v1_noNoGroup_metadata.csv             — obs metadata (for 06_h5ad_to_seurat.R)
  scib_lum_v1_30pc3000g/scib_results_*.csv       — benchmark results

Dependencies: scanpy, scvi-tools, bbknn, scanorama, scib-metrics
"""

import os
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
import bbknn
import scanorama
import scvi
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scib_metrics.benchmark import Benchmarker

sc.settings.verbosity = 3
sc.logging.print_header()

# =============================================================================
# Paths
# =============================================================================
BASE_DIR       = "/mctp/share/users/hanbyul/projects/human_prostate_sc/cellranger_7.1/17datasets/python_final/032025_scib"
PLOT_DIR       = os.path.join(BASE_DIR, "scib_lum_v1_30pc3000g")
CHECKPOINT_DIR = os.path.join(BASE_DIR, "checkpoints_lum_v1_30pc3000g")
H5AD_IN        = os.path.join(BASE_DIR, "seu_lum_v1_noNoGroup.h5ad")         # from 03_batch_correction_luminal.R
H5AD_OUT       = os.path.join(BASE_DIR, "seu_lum_v1_noNoGroup_scib_30pc3000g.h5ad")
OBS_CSV        = os.path.join(BASE_DIR, "seu_lum_v1_noNoGroup_metadata.csv")  # full metadata from step 02

os.makedirs(PLOT_DIR,       exist_ok=True)
os.makedirs(CHECKPOINT_DIR, exist_ok=True)
sc.settings.figdir = PLOT_DIR

BATCH_KEY = "Paper_ID_short"
LABEL_KEY = "etsgroup"

N_HVG       = 3000
N_PCS       = 30
N_NEIGHBORS = 30

# =============================================================================
# Checkpoint helpers
# =============================================================================
def save_checkpoint(adata, tag):
    path = os.path.join(CHECKPOINT_DIR, f"ckpt_{tag}.h5ad")
    print(f"[checkpoint] saving → {path}")
    adata.write_h5ad(path)

def load_checkpoint(tag):
    path = os.path.join(CHECKPOINT_DIR, f"ckpt_{tag}.h5ad")
    if os.path.exists(path):
        print(f"[checkpoint] loading ← {path}")
        return sc.read_h5ad(path)
    return None

# =============================================================================
# Color palettes
# =============================================================================
paper_col = {
    "P01":"#009CD2","P02":"#e6194B","P03":"#3cb44b","P05":"#ffe119","P06":"#4363d8",
    "P07":"#f58231","P08":"#911eb4","P09":"#46f0f0","P10":"#2b1700","P11":"#f032e6",
    "P12":"#bcf60c","P14":"#fabebe","P15":"#008080","P16":"#e6beff","P17":"#9A6324",
    "P19":"#FFFF33","P21":"#999999",
}
etsgroup_col = {
    "N":"#71A45E","Benign":"#87CEFA","BPH":"#0047ab","nonETS":"#F781BF",
    "ERG":"#984EA3","ETV":"#FF7F00","AR":"#E41A1C","NEPC":"#A65628","DN":"#4DAF4A",
}

def get_palette(adata, key, col_dict):
    cats = adata.obs[key].cat.categories.tolist()
    return [col_dict.get(c, "#AAAAAA") for c in cats]

# =============================================================================
# Step 1: Load h5ad exported by 03_batch_correction_luminal.R
# The h5ad was saved by SeuratDisk (DietSeurat → SaveH5Seurat → Convert).
# SeuratDisk puts log-normalized data in raw slot and scale.data in X.
# We restore lognorm from raw, back-compute pseudo-counts, and re-attach
# the full obs metadata from the CSV.
# =============================================================================
print("\n" + "="*60)
print("Step 1: Load Seurat h5ad (output of 03_batch_correction_luminal.R)")
print("="*60)

ckpt = load_checkpoint("after_load")
if ckpt is not None:
    adata = ckpt
else:
    adata = sc.read_h5ad(H5AD_IN)
    print(adata)
    print("obsm keys on load:", list(adata.obsm.keys()))

    # Preserve Seurat reductions before raw.to_adata() wipes obsm
    obsm_backup = dict(adata.obsm)

    # Restore log-normalized matrix and back-compute pseudo-counts for scVI
    adata = adata.raw.to_adata()
    adata.layers["lognorm"] = adata.X.copy()
    raw_dense  = adata.X.toarray() if sp.issparse(adata.X) else np.array(adata.X)
    adata.layers["counts"] = sp.csr_matrix(np.expm1(raw_dense).round().astype(np.float32))

    # Restore Seurat embeddings (CCA, RPCA, Harmony, UMAP variants)
    for k, v in obsm_backup.items():
        adata.obsm[k] = v

    # Re-attach full obs metadata from CSV
    # (metadata was separated from the stripped h5ad in 02_data_merge_and_filter.py)
    cell_df = pd.read_csv(OBS_CSV, index_col=0).reindex(adata.obs_names)
    new_cols = cell_df.columns.difference(adata.obs.columns)
    if len(new_cols) > 0:
        adata.obs = pd.concat([adata.obs, cell_df[new_cols]], axis=1)

    for col in adata.obs.select_dtypes(include="object").columns:
        adata.obs[col] = adata.obs[col].astype(str)
    for key in [BATCH_KEY, LABEL_KEY]:
        if key in adata.obs.columns:
            adata.obs[key] = adata.obs[key].fillna("Unknown").astype("category")
        else:
            raise KeyError(f"Column {key!r} not found. Check obs_csv path.")

    print(adata)
    save_checkpoint(adata, "after_load")

# =============================================================================
# Step 2: Rename Seurat embedding keys to standard Python convention
# Seurat exports reduction names with dots (e.g., "X_umap.cca");
# we rename to underscores for clean obsm key names.
# =============================================================================
print("\n" + "="*60)
print("Step 2: Rename Seurat embedding keys")
print("="*60)

rename_map = {
    # Seurat name (as in h5ad)        : Python standard name
    "X_pca"               : "X_pca_unintegrated",    # unintegrated PCA baseline
    "X_umap.unintegrated" : "X_umap_unintegrated",
    "X_integrated.cca"    : "X_cca",                  # CCA corrected embedding
    "X_umap.cca"          : "X_umap_cca",
    "X_integrated.rpca"   : "X_rpca",                 # RPCA corrected embedding
    "X_umap.rpca"         : "X_umap_rpca",
    "X_harmony"           : "X_harmony",
    "X_umap.harmony"      : "X_umap_harmony",
}
for old, new in rename_map.items():
    if old in adata.obsm:
        adata.obsm[new] = adata.obsm[old]
        if old != new:
            del adata.obsm[old]
        print(f"  {old} → {new}")

print("obsm keys after rename:", list(adata.obsm.keys()))

# =============================================================================
# Step 3: HVG selection + Python PCA
# (Used as input for Scanorama and BBKNN; also serves as the Python-computed
#  unintegrated baseline for Scanorama and BBKNN initialization)
# =============================================================================
print("\n" + "="*60)
print(f"Step 3: HVG ({N_HVG}) + Python PCA ({N_PCS} PCs)")
print("="*60)

ckpt = load_checkpoint("after_pca")
if ckpt is not None:
    adata = ckpt
    adata_hvg = None
else:
    sc.pp.highly_variable_genes(
        adata, n_top_genes=N_HVG, batch_key=BATCH_KEY,
        flavor="seurat_v3", layer="counts"
    )
    adata_hvg = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata_hvg)
    sc.pp.pca(adata_hvg, n_comps=N_PCS)
    adata.obsm["X_pca_python"] = adata_hvg.obsm["X_pca"].copy()
    print("Python PCA:", adata.obsm["X_pca_python"].shape)
    save_checkpoint(adata, "after_pca")

# =============================================================================
# Step 4: Scanorama (dimred=30, matching N_PCS for consistent output dimensions)
# =============================================================================
print("\n" + "="*60)
print("Step 4: Scanorama (dimred=30)")
print("="*60)

ckpt = load_checkpoint("after_scanorama")
if ckpt is not None:
    adata = ckpt
else:
    if adata_hvg is None:
        sc.pp.highly_variable_genes(adata, n_top_genes=N_HVG, batch_key=BATCH_KEY,
                                    flavor="seurat_v3", layer="counts")
        adata_hvg = adata[:, adata.var["highly_variable"]].copy()
        sc.pp.scale(adata_hvg); sc.pp.pca(adata_hvg, n_comps=N_PCS)
        adata.obsm["X_pca_python"] = adata_hvg.obsm["X_pca"].copy()

    batches = sorted(adata_hvg.obs[BATCH_KEY].unique())
    adatas_list = [adata_hvg[adata_hvg.obs[BATCH_KEY] == b].copy() for b in batches]
    scanorama.integrate_scanpy(adatas_list, dimred=N_PCS)
    emb_df = pd.concat([pd.DataFrame(a.obsm["X_scanorama"], index=a.obs_names)
                         for a in adatas_list])
    adata.obsm["X_scanorama"] = emb_df.loc[adata.obs_names].values
    print("Scanorama:", adata.obsm["X_scanorama"].shape)
    save_checkpoint(adata, "after_scanorama")

# =============================================================================
# Step 5: BBKNN (graph-based; UMAP used as proxy for scib)
# =============================================================================
print("\n" + "="*60)
print("Step 5: BBKNN")
print("="*60)

ckpt = load_checkpoint("after_bbknn")
if ckpt is not None:
    adata = ckpt
else:
    adata_bbknn = adata.copy()
    adata_bbknn.obsm["X_pca"] = adata.obsm["X_pca_python"].copy()
    bbknn.bbknn(adata_bbknn, batch_key=BATCH_KEY, n_pcs=N_PCS)
    sc.tl.umap(adata_bbknn)
    adata.obsm["X_umap_bbknn"] = adata_bbknn.obsm["X_umap"].copy()
    adata.obsp["bbknn_connectivities"] = adata_bbknn.obsp["connectivities"].copy()
    save_checkpoint(adata, "after_bbknn")

# =============================================================================
# Step 6: scVI + scANVI (3,000 HVG subset, pseudo-counts as input)
# Key update vs earlier notebook: scVI/scANVI are trained on the same 3,000 HVG
# subset used by Scanorama and BBKNN (earlier notebook used full count matrix).
# This ensures all methods operate on the same gene set for fair comparison.
# scANVI uses ETS subgroup labels (etsgroup) for label-guided integration.
# =============================================================================
print("\n" + "="*60)
print("Step 6: scVI + scANVI (3000 HVG subset)")
print("="*60)

ckpt = load_checkpoint("after_scvi_scanvi")
if ckpt is not None:
    adata = ckpt
else:
    if "highly_variable" not in adata.var.columns:
        sc.pp.highly_variable_genes(adata, n_top_genes=N_HVG, batch_key=BATCH_KEY,
                                    flavor="seurat_v3", layer="counts")

    adata_scvi = adata[:, adata.var["highly_variable"]].copy()
    adata_scvi.X = adata_scvi.layers["counts"].copy()   # pseudo-counts as input

    scvi.model.SCVI.setup_anndata(adata_scvi, batch_key=BATCH_KEY)
    model_scvi = scvi.model.SCVI(adata_scvi, n_latent=30, n_layers=2)
    model_scvi.train(max_epochs=400, early_stopping=True,
                     datasplitter_kwargs={"num_workers": 0})
    adata.obsm["X_scvi"] = model_scvi.get_latent_representation()

    adata_scvi.obs["celltype_scanvi"] = adata.obs[LABEL_KEY].astype(str).values
    model_scanvi = scvi.model.SCANVI.from_scvi_model(
        model_scvi, adata=adata_scvi,
        labels_key="celltype_scanvi", unlabeled_category="Unknown"
    )
    model_scanvi.train(max_epochs=20, n_samples_per_label=100,
                       datasplitter_kwargs={"num_workers": 0})
    adata.obsm["X_scanvi"] = model_scanvi.get_latent_representation()
    print("scVI:", adata.obsm["X_scvi"].shape, "  scANVI:", adata.obsm["X_scanvi"].shape)
    save_checkpoint(adata, "after_scvi_scanvi")

# =============================================================================
# Step 7: UMAP embeddings for Python methods
# (CCA, RPCA, Harmony, unintegrated UMAPs already in obsm from R step)
# =============================================================================
print("\n" + "="*60)
print("Step 7: UMAP embeddings")
print("="*60)

ckpt = load_checkpoint("after_umaps")
if ckpt is not None:
    adata = ckpt
else:
    def compute_umap(adata, embed_key, umap_key, n_neighbors=N_NEIGHBORS):
        tmp = adata.copy()
        tmp.obsm["X_emb"] = adata.obsm[embed_key]
        sc.pp.neighbors(tmp, use_rep="X_emb", n_neighbors=n_neighbors)
        sc.tl.umap(tmp)
        adata.obsm[umap_key] = tmp.obsm["X_umap"].copy()
        print(f"  {umap_key}: done")

    compute_umap(adata, "X_scanorama", "X_umap_scanorama")
    compute_umap(adata, "X_scvi",      "X_umap_scvi")
    compute_umap(adata, "X_scanvi",    "X_umap_scanvi")
    # X_umap_unintegrated, X_umap_cca, X_umap_rpca, X_umap_harmony → from R
    # X_umap_bbknn → from Step 5
    save_checkpoint(adata, "after_umaps")

# =============================================================================
# Step 8: UMAP visualization
# =============================================================================
print("\n" + "="*60)
print("Step 8: UMAP plots")
print("="*60)

adata.obs["Paper_ID_short2"] = (
    adata.obs[BATCH_KEY].astype(str).str.extract(r"(P\d+)")[0].astype("category")
)
methods_umap = [
    ("X_umap_unintegrated","Unintegrated (PCA)"), ("X_umap_cca","CCA (Seurat)"),
    ("X_umap_rpca","RPCA (Seurat)"),             ("X_umap_harmony","Harmony (Seurat)"),
    ("X_umap_scanorama","Scanorama"),             ("X_umap_bbknn","BBKNN"),
    ("X_umap_scvi","scVI"),                       ("X_umap_scanvi","scANVI"),
]
methods_umap = [(k, t) for k, t in methods_umap if k in adata.obsm]

fig, axes = plt.subplots(len(methods_umap), 2, figsize=(14, 5 * len(methods_umap)))
if len(methods_umap) == 1:
    axes = axes[np.newaxis, :]
for row, (umap_key, title) in enumerate(methods_umap):
    sc.pl.embedding(adata, basis=umap_key, color="Paper_ID_short2",
                    palette=get_palette(adata, "Paper_ID_short2", paper_col),
                    ax=axes[row, 0], show=False, frameon=False,
                    title=f"{title} — Batch", s=3)
    sc.pl.embedding(adata, basis=umap_key, color=LABEL_KEY,
                    palette=get_palette(adata, LABEL_KEY, etsgroup_col),
                    ax=axes[row, 1], show=False, frameon=False,
                    title=f"{title} — ETS subgroup", s=3)
plt.tight_layout()
fig.savefig(os.path.join(PLOT_DIR, "umap_all_methods.pdf"), bbox_inches="tight")
fig.savefig(os.path.join(PLOT_DIR, "umap_all_methods.png"), dpi=100, bbox_inches="tight")
plt.close("all")

# =============================================================================
# Step 9: scib-metrics benchmarking (Figure 5b, d; Supplementary Tables 5, 6)
# Label key: etsgroup (ETS molecular subgroup) for bio conservation
# =============================================================================
print("\n" + "="*60)
print("Step 9: scib-metrics benchmark (Figure 5b, d)")
print("="*60)

embed_keys = [
    "X_pca_unintegrated",   # unintegrated baseline (Seurat PCA)
    "X_cca",                # CCA (Seurat)
    "X_rpca",               # RPCA (Seurat)
    "X_harmony",            # Harmony (Seurat)
    "X_scanorama",          # Scanorama (Python)
    "X_umap_bbknn",         # BBKNN proxy (Python)
    "X_scvi",               # scVI (Python)
    "X_scanvi",             # scANVI (Python) ← selected
]
embed_keys = [k for k in embed_keys if k in adata.obsm]
print("Benchmarking:", embed_keys)

LABEL_MAP = {
    "X_pca_unintegrated":"Unintegrated","X_cca":"CCA","X_rpca":"RPCA",
    "X_harmony":"Harmony","X_scanorama":"Scanorama","X_umap_bbknn":"BBKNN†",
    "X_scvi":"scVI","X_scanvi":"scANVI",
}

bm = Benchmarker(adata, batch_key=BATCH_KEY, label_key=LABEL_KEY,
                 embedding_obsm_keys=embed_keys, n_jobs=1)
bm.benchmark()

def rename_methods(df):
    df = df.copy()
    df.index = [LABEL_MAP.get(i, i) for i in df.index]
    return df

df_raw    = rename_methods(bm.get_results(min_max_scale=False))
df_scaled = rename_methods(bm.get_results(min_max_scale=True))

df_raw.to_csv(   os.path.join(PLOT_DIR, "scib_results_raw.csv"))
df_scaled.to_csv(os.path.join(PLOT_DIR, "scib_results_minmax.csv"))
print("Saved: scib_results_raw.csv  (Supplementary Table 5)")
print("Saved: scib_results_minmax.csv  (Supplementary Table 6)")
print(df_scaled)

for scale, tag in [(False, "_raw"), (True, "_minmax")]:
    fig_t = bm.plot_results_table(min_max_scale=scale, show=False)
    fig_t.figure.savefig(os.path.join(PLOT_DIR, f"scib_results_table{tag}.pdf"), bbox_inches="tight")
    plt.close("all")

def make_scatter(df, title, out_stem):
    batch_col, bio_col = "Batch correction", "Bio conservation"
    df_p = df[df.index != "Metric Type"].copy()
    for c in [batch_col, bio_col]:
        df_p[c] = pd.to_numeric(df_p[c], errors="coerce")
    df_p = df_p.dropna(subset=[batch_col, bio_col])
    fig, ax = plt.subplots(figsize=(6, 5))
    for i, (method, row) in enumerate(df_p.iterrows()):
        is_bbknn = "BBKNN" in method
        ax.scatter(row[batch_col], row[bio_col],
                   s=160 if is_bbknn else 100, marker="*" if is_bbknn else "o",
                   color=plt.cm.tab10.colors[i % 10], zorder=3,
                   edgecolors="k" if is_bbknn else "none", linewidths=0.5)
        ax.annotate(method, (row[batch_col], row[bio_col]),
                    textcoords="offset points", xytext=(6, 3), fontsize=8)
    ax.set_xlabel("Batch correction", fontsize=12); ax.set_ylabel("Bio conservation", fontsize=12)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.axhline(0.5, color="gray", lw=0.8, ls="--"); ax.axvline(0.5, color="gray", lw=0.8, ls="--")
    ax.set_title(title, fontsize=10)
    plt.tight_layout()
    fig.savefig(os.path.join(PLOT_DIR, f"{out_stem}.pdf"), bbox_inches="tight")
    plt.close("all")

make_scatter(df_scaled, "Batch correction vs Bio conservation (min-max scaled)\nLuminal cells", "scib_scatter_minmax")

# =============================================================================
# Step 10: Save final h5ad
# =============================================================================
print("\n" + "="*60)
print("Step 10: Save final h5ad")
print("="*60)

adata.X = adata.layers["lognorm"].copy()
adata.write_h5ad(H5AD_OUT)
print(f"Saved: {H5AD_OUT}")
print(f"→ Next: run 06_h5ad_to_seurat.R to convert to Seurat RDS for plotting")
