"""
04_integration_all_cells.py
============================
Step 4 of the pipeline: batch correction and benchmarking for all cell types.

Cho et al. — A harmonized single-cell RNA-seq atlas of the human prostate
Scientific Data (in submission)

Workflow position:
  02_data_merge_and_filter.py → [THIS FILE] → 06_h5ad_to_seurat.R

Input:
  merged_a17_hit1.h5ad   — created by 02_data_merge_and_filter.py
  (all annotated cells, low-quality cells removed, MT genes removed,
   raw counts in layers["counts"], log-normalized in layers["lognorm"])

Parameters (manuscript: all cell types):
  - 2,000 highly variable genes (batch-aware, Seurat v3 flavor, raw counts)
  - 50 principal components
  - Batch key : Paper_ID_short (study short code, e.g. P01)
  - Label key : final_hit1_group_detailed (annotated cell type)
  - kNN neighbors for UMAP: 30

Integration methods (6 total, all benchmarked):
  1. Unintegrated (PCA)  — baseline
  2. Harmony             — iterative cluster-based correction on PCA embedding
  3. Scanorama           — panoramic stitching in PCA space (dimred=50)
  4. BBKNN               — graph-only; UMAP used as proxy for benchmarking †
  5. scVI                — n_latent=30, n_layers=2, max_epochs=400, early stopping
  6. scANVI              — n_latent=30, max_epochs=20, n_samples_per_label=100 ← selected

† BBKNN corrects the neighbor graph, not a low-dimensional embedding.
  The UMAP derived from the BBKNN graph is used as a proxy for scib.
  PCR comparison metric (requires corrected high-dim embedding) is excluded
  from BBKNN's aggregate score.

Benchmarking: scib-metrics (Luecken et al. 2022, Nat Methods)
  Total = 0.4 × Batch correction + 0.6 × Bio conservation
  Results → scib_results_raw.csv / scib_results_minmax.csv (Supplementary Tables 3, 4)

Output:
  merged_a17_hit1_scib_50pc2000g.h5ad       — all embeddings, X = lognorm
  merged_a17_hit1_metadata.csv              — obs metadata (for 06_h5ad_to_seurat.R)
  scib_a17_hit1_50pc2000g/scib_results_*.csv — benchmark results

Dependencies: scanpy, scvi-tools, bbknn, scanorama, harmonypy, scib-metrics
"""

import os
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
import bbknn
import scanorama
import scvi
import harmonypy as hm
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scib_metrics.benchmark import Benchmarker

sc.settings.verbosity = 3
sc.logging.print_header()

# =============================================================================
# Paths
# =============================================================================
BASE_DIR       = "/mctp/share/users/hanbyul/projects/human_prostate_sc/cellranger_7.1/17datasets/python_final"
PLOT_DIR       = os.path.join(BASE_DIR, "scib_a17_hit1_50pc2000g")
CHECKPOINT_DIR = os.path.join(BASE_DIR, "checkpoints_a17_hit1_50pc2000g")
H5AD_IN        = os.path.join(BASE_DIR, "merged_a17_hit1.h5ad")
H5AD_OUT       = os.path.join(BASE_DIR, "merged_a17_hit1_scib_50pc2000g.h5ad")
META_CSV_OUT   = os.path.join(BASE_DIR, "merged_a17_hit1_metadata.csv")   # for 06_h5ad_to_seurat.R

os.makedirs(PLOT_DIR,       exist_ok=True)
os.makedirs(CHECKPOINT_DIR, exist_ok=True)
sc.settings.figdir = PLOT_DIR

BATCH_KEY    = "Paper_ID_short"
LABEL_KEY    = "final_hit1_group_detailed"
ETSGROUP_KEY = "etsgroup"

N_HVG       = 2000
N_PCS       = 50
N_NEIGHBORS = 30

# =============================================================================
# Checkpoint helpers (resume from any step after interruption)
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
group_detailed_col = {
    "Luminal":"#E31A1C","Club":"#E2ADB2","Basal":"#1976D2","Hillock":"#87CEFA",
    "Ionocyte":"#4A1486","T cell":"#3B7A57","B cell":"#4DAF4A","Plasma":"#008B8B",
    "NK cell":"#00F0A8","Granulocytes":"#8FBC8F","Macrophage":"#00AB66",
    "Monocyte":"#CEFF00","Dendritic":"#5D4037","Endothelial cell":"#FFD800",
    "SMC":"#FD8D3C","Fibroblast":"#F0E68C","Schwann":"#696969",
}

def get_palette(adata, key, col_dict):
    cats = adata.obs[key].cat.categories.tolist()
    return [col_dict.get(c, "#AAAAAA") for c in cats]

# =============================================================================
# Step 1: Load and prepare
# (merged_a17_hit1.h5ad already has counts/lognorm layers from step 02)
# =============================================================================
print("\n" + "="*60)
print("Step 1: Load merged_a17_hit1.h5ad")
print("="*60)

ckpt = load_checkpoint("after_load")
if ckpt is not None:
    adata = ckpt
else:
    adata = sc.read_h5ad(H5AD_IN)
    print(adata)
    print(f"Cells: {adata.n_obs:,}  Genes: {adata.n_vars:,}")

    # Confirm layers
    if "counts" not in adata.layers:
        # If loaded from merged_a17_hit1.h5ad which has X = raw counts
        adata.layers["counts"] = adata.X.copy()
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        adata.layers["lognorm"] = adata.X.copy()
    else:
        adata.X = adata.layers["lognorm"].copy()

    for col in adata.obs.select_dtypes(include="object").columns:
        adata.obs[col] = adata.obs[col].astype(str)
    for key in [BATCH_KEY, LABEL_KEY, ETSGROUP_KEY]:
        if key in adata.obs.columns:
            adata.obs[key] = adata.obs[key].fillna("Unknown").astype("category")

    save_checkpoint(adata, "after_load")

# Export obs metadata CSV for 06_h5ad_to_seurat.R
adata.obs.to_csv(META_CSV_OUT)
print(f"Saved obs metadata: {META_CSV_OUT}")

# =============================================================================
# Step 2: HVG selection + PCA (unintegrated baseline)
# =============================================================================
print("\n" + "="*60)
print(f"Step 2: HVG ({N_HVG}) + PCA ({N_PCS} PCs) — unintegrated baseline")
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
    adata.obsm["X_pca_unintegrated"] = adata_hvg.obsm["X_pca"].copy()
    print("Unintegrated PCA:", adata.obsm["X_pca_unintegrated"].shape)
    save_checkpoint(adata, "after_pca")

# =============================================================================
# Step 3: Harmony
# =============================================================================
print("\n" + "="*60)
print("Step 3: Harmony")
print("="*60)

ckpt = load_checkpoint("after_harmony")
if ckpt is not None:
    adata = ckpt
else:
    ho = hm.run_harmony(adata.obsm["X_pca_unintegrated"], adata.obs, BATCH_KEY)
    adata.obsm["X_harmony"] = ho.Z_corr.T
    print("Harmony:", adata.obsm["X_harmony"].shape)
    save_checkpoint(adata, "after_harmony")

# =============================================================================
# Step 4: Scanorama
# =============================================================================
print("\n" + "="*60)
print("Step 4: Scanorama")
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
    batches = sorted(adata_hvg.obs[BATCH_KEY].cat.categories)
    adatas_list = [adata_hvg[adata_hvg.obs[BATCH_KEY] == b].copy() for b in batches]
    scanorama.integrate_scanpy(adatas_list, dimred=N_PCS)
    emb_df = pd.concat([pd.DataFrame(a.obsm["X_scanorama"], index=a.obs_names)
                         for a in adatas_list])
    adata.obsm["X_scanorama"] = emb_df.loc[adata.obs_names].values
    print("Scanorama:", adata.obsm["X_scanorama"].shape)
    save_checkpoint(adata, "after_scanorama")

# =============================================================================
# Step 5: BBKNN
# BBKNN corrects the neighbor graph (not an embedding).
# UMAP from the corrected graph is used as a proxy for scib benchmarking.
# PCR comparison is excluded from BBKNN's aggregate score (see manuscript methods).
# =============================================================================
print("\n" + "="*60)
print("Step 5: BBKNN (graph-based; UMAP used as proxy for scib)")
print("="*60)

ckpt = load_checkpoint("after_bbknn")
if ckpt is not None:
    adata = ckpt
else:
    adata_bbknn = adata.copy()
    adata_bbknn.obsm["X_pca"] = adata.obsm["X_pca_unintegrated"].copy()
    bbknn.bbknn(adata_bbknn, batch_key=BATCH_KEY, n_pcs=N_PCS)
    sc.tl.umap(adata_bbknn)
    adata.obsm["X_umap_bbknn"] = adata_bbknn.obsm["X_umap"].copy()
    adata.obsp["bbknn_connectivities"] = adata_bbknn.obsp["connectivities"].copy()
    save_checkpoint(adata, "after_bbknn")

# =============================================================================
# Step 6: scVI + scANVI (trained on 2,000 HVG subset, raw counts)
# scANVI is label-guided (uses cell type annotations from R step 01).
# scANVI was selected as the primary integration method based on benchmarking.
# =============================================================================
print("\n" + "="*60)
print("Step 6: scVI + scANVI")
print("="*60)

ckpt = load_checkpoint("after_scvi_scanvi")
if ckpt is not None:
    adata = ckpt
else:
    if "highly_variable" not in adata.var.columns:
        sc.pp.highly_variable_genes(adata, n_top_genes=N_HVG, batch_key=BATCH_KEY,
                                    flavor="seurat_v3", layer="counts")
    adata_scvi = adata[:, adata.var["highly_variable"]].copy()
    adata_scvi.X = adata_scvi.layers["counts"].copy()

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
# Step 7: UMAP embeddings
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

    compute_umap(adata, "X_pca_unintegrated", "X_umap_unintegrated")
    compute_umap(adata, "X_harmony",          "X_umap_harmony")
    compute_umap(adata, "X_scanorama",        "X_umap_scanorama")
    compute_umap(adata, "X_scvi",             "X_umap_scvi")
    compute_umap(adata, "X_scanvi",           "X_umap_scanvi")
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
    ("X_umap_unintegrated","Unintegrated (PCA)"), ("X_umap_harmony","Harmony"),
    ("X_umap_scanorama","Scanorama"), ("X_umap_bbknn","BBKNN"),
    ("X_umap_scvi","scVI"), ("X_umap_scanvi","scANVI"),
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
                    palette=get_palette(adata, LABEL_KEY, group_detailed_col),
                    ax=axes[row, 1], show=False, frameon=False,
                    title=f"{title} — Cell type", s=3)
plt.tight_layout()
fig.savefig(os.path.join(PLOT_DIR, "umap_all_methods.pdf"), bbox_inches="tight")
fig.savefig(os.path.join(PLOT_DIR, "umap_all_methods.png"), dpi=100, bbox_inches="tight")
plt.close("all")

# =============================================================================
# Step 9: scib-metrics benchmarking (Figure 5a, c; Supplementary Tables 3, 4)
# =============================================================================
print("\n" + "="*60)
print("Step 9: scib-metrics benchmark (Figure 5a, c)")
print("="*60)

embed_keys = ["X_pca_unintegrated","X_harmony","X_scanorama",
              "X_umap_bbknn","X_scvi","X_scanvi"]
embed_keys = [k for k in embed_keys if k in adata.obsm]

LABEL_MAP = {
    "X_pca_unintegrated":"Unintegrated","X_harmony":"Harmony",
    "X_scanorama":"Scanorama","X_umap_bbknn":"BBKNN†",
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
print("Saved: scib_results_raw.csv  (Supplementary Table 3)")
print("Saved: scib_results_minmax.csv  (Supplementary Table 4)")
print(df_scaled)

for scale, tag in [(False, "_raw"), (True, "_minmax")]:
    fig_t = bm.plot_results_table(min_max_scale=scale, show=False)
    fig_t.figure.savefig(os.path.join(PLOT_DIR, f"scib_results_table{tag}.pdf"), bbox_inches="tight")
    plt.close("all")

# Scatter plot (Figure 5a)
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
    ax.set_xlabel("Batch correction", fontsize=12)
    ax.set_ylabel("Bio conservation", fontsize=12)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.axhline(0.5, color="gray", lw=0.8, ls="--"); ax.axvline(0.5, color="gray", lw=0.8, ls="--")
    ax.set_title(title, fontsize=10)
    fig.text(0.01, -0.02, "† BBKNN: proxy score via X_umap_bbknn (2D UMAP)", fontsize=7, color="gray")
    plt.tight_layout()
    fig.savefig(os.path.join(PLOT_DIR, f"{out_stem}.pdf"), bbox_inches="tight")
    plt.close("all")

make_scatter(df_scaled, "Batch correction vs Bio conservation (min-max scaled)\nAll cell types", "scib_scatter_minmax")

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
