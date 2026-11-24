# Code release header: re-scores the stored luminal embeddings with the scib-metrics Benchmarker on the released labels. No integration is re-run.
# Input: one stored luminal run folder (input object, embeddings, run_settings.json) and the cell label table (--label-table, md5 checked).
# Output: scib_results_*.csv and rescore_settings.json.
# Cells whose label is in --unlabelled-values (Unassigned for the luminal benchmark) are not scored. Seed, batch key, threads and package versions must equal the run's record.
# It loads umap_sigma_fix.py from the run's script folder, with the same code as 05_integration/umap_sigma_fix.py. Methods: Benchmarking of integration methods.
"""R1_rescore_scib.py - re-score the stored Luminal-object embeddings (final luminal integration run) with the scIB Benchmarker,
changing only the label column. No integration is re-run.

What is the SAME as the run (the run's copy of 03_run_python_integration.py, md5 a1a0c3f97ef9529c71f8939555e98f36):
  * thread guard (visible cores = --threads = every thread variable), imports in the same order, torch threads,
    the umap sigma fix (the run's umap_sigma_fix.py, loaded from the run's script folder, guard() must pass);
  * step 1 (load input.h5ad of the run, integer counts check, normalize_total 1e4, log1p, batch and label columns
    as str -> category with fillna 'Unknown'), the R embedding load (r_embeddings.h5, same reorder), hvg_var.csv into
    adata.var as 03 does on a resume;
  * the stored Python embeddings (python/embeddings/*.npy), as 03 loads them on a resume;
  * step 9 code: seed_all(seed of the run); subset to the scored cells (adata[_keep].copy(), unused categories
    removed); Benchmarker(bdata, batch_key, label_key, embedding_obsm_keys, n_jobs=1) with package defaults;
    get_results(min_max_scale=True/False) and _results written as the run did.
What changes: only the label column and therefore which cells are 'labeled':
  --labels old   : the run's own label column (luminal_scanvi_label), scored = label != 'Unknown' (reproduces
                   the stored tables; used for the small test);
  --labels final : release_subgroup_value of the label table (luminal_cells_v7_base_v25.csv.gz), joined by
                   cell_id; scored = final label not in --unlabelled-values (default 'Unassigned,Not assessed').
The run's settings (versions, sigma-fix md5, seed, batch key, threads) are read from python/run_settings.json and the
script stops if this process differs.
"""
import argparse, hashlib, json, os, platform, random, re, resource, sys, time

p = argparse.ArgumentParser()
p.add_argument('--run-dir', required=True, help='reintegration_v18/runs/<run>/luminal_seedN')
p.add_argument('--outdir', required=True)
p.add_argument('--labels', choices=['old', 'final'], required=True)
p.add_argument('--label-table', default=None)
p.add_argument('--label-table-md5', default=None)
p.add_argument('--final-column', default='release_subgroup_value')
p.add_argument('--unlabelled-values', default='Unassigned,Not assessed')
p.add_argument('--embeddings', default=None, help='comma list of obsm keys (default: the run\'s 8, in its order)')
p.add_argument('--threads', type=int, required=True)
a = p.parse_args()

V18 = '<PROJECT_DIR>/reintegration_v18/scripts'
RUN = os.path.abspath(a.run_dir)
ST = json.load(open(os.path.join(RUN, 'python', 'run_settings.json')))
S0, V0 = ST['settings'], ST['versions']
assert S0['mode'] == 'luminal' and S0['scib_cells'] == 'labelled', (S0['mode'], S0['scib_cells'])

# 03's thread guard
_aff = len(os.sched_getaffinity(0))
_tv = {k: os.environ.get(k) for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMBA_NUM_THREADS')}
_bad = {k: v for k, v in _tv.items() if v != str(a.threads)}
if _aff != a.threads or _bad:
    sys.exit(f'--threads {a.threads} but {_aff} visible cores and thread variables {_tv}')
if a.threads != S0['threads']:
    sys.exit(f'--threads {a.threads} but the run used {S0["threads"]}')

MODE = 'luminal'
OLD_LABEL_KEY = S0['label_key']           # luminal_scanvi_label
BATCH_KEY = S0['batch_key']               # Paper_ID_short
SEED = int(S0['seed'])
LABEL_KEY = OLD_LABEL_KEY if a.labels == 'old' else 'final_subgroup_label'

# imports in 03's order
import numpy as np
import pandas as pd
import scipy.sparse as sp
import h5py
import anndata as ad
import scanpy as sc
import torch
torch.set_num_threads(a.threads)
import bbknn
import scanorama
import scvi
import harmonypy as hm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scib_metrics.benchmark import Benchmarker
sys.dont_write_bytecode = True; sys.path.insert(0, V18)
import umap_sigma_fix
del sys.path[0]
SIGMA_FIX = umap_sigma_fix.guard()
if SIGMA_FIX.get('module_md5') != S0['umap_sigma_fix_md5']:
    sys.exit(f'sigma fix md5 {SIGMA_FIX.get("module_md5")} != run {S0["umap_sigma_fix_md5"]}')

import importlib.metadata as im
VERS = {k: im.version(k) for k in ['scanpy', 'anndata', 'scvi-tools', 'scib-metrics', 'harmonypy', 'scanorama', 'bbknn',
                                   'torch', 'lightning', 'jax', 'jaxlib', 'numpy', 'scipy', 'pandas', 'scikit-learn',
                                   'umap-learn', 'pynndescent', 'numba', 'llvmlite', 'scikit-misc', 'h5py']}
VERS.update(python=platform.python_version(), executable=sys.executable)
_vd = {k: (V0.get(k), v) for k, v in VERS.items() if V0.get(k) != v}
if _vd:
    sys.exit(f'package versions differ from the run: {_vd}')

OUT = os.path.abspath(a.outdir)
os.makedirs(OUT, exist_ok=True)
T0 = time.time()
def log(*x):
    print(time.strftime('%Y-%m-%d %H:%M:%S'), f'[rescore +{(time.time() - T0) / 60:.1f} min]', *x, flush=True)
def seed_all(s=SEED):
    random.seed(s); np.random.seed(s); torch.manual_seed(s)
def fsync_replace(tmp, final):
    with open(tmp, 'rb+') as fh:
        os.fsync(fh.fileno())
    if os.path.getsize(tmp) == 0:
        sys.exit(f'{tmp} is empty after writing')
    os.replace(tmp, final)
def write_csv(df, final):
    tmp = f'{final}.tmp{os.getpid()}'; df.to_csv(tmp); fsync_replace(tmp, final)
def write_json(obj, final):
    tmp = f'{final}.tmp{os.getpid()}'
    with open(tmp, 'w') as fh:
        json.dump(obj, fh, indent=1)
    fsync_replace(tmp, final)
def md5(path):
    h = hashlib.md5()
    with open(path, 'rb') as fh:
        for b in iter(lambda: fh.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()

EMB = os.path.join(RUN, 'python', 'embeddings')
INPUT = S0['input']
assert os.path.abspath(INPUT) == os.path.join(RUN, 'input.h5ad'), INPUT
REC = dict(run_dir=RUN, run_id=S0['run_id'], seed=SEED, labels=a.labels, label_key=LABEL_KEY, batch_key=BATCH_KEY,
           threads=a.threads, cpu_affinity=sorted(os.sched_getaffinity(0)), thread_env=_tv, host=platform.node(),
           cpu_model=next((l.split(':', 1)[1].strip() for l in open('/proc/cpuinfo') if l.startswith('model name')), None),
           umap_sigma_fix=SIGMA_FIX, versions=VERS, started=time.strftime('%Y-%m-%d %H:%M:%S'),
           script=os.path.abspath(__file__), script_md5=md5(os.path.abspath(__file__)))

# ── step 1 as 03 ──
adata = ad.read_h5ad(INPUT)
assert sp.issparse(adata.X)
d = adata.X.data
assert np.all(d == np.floor(d)), 'X must hold integer raw counts'
adata.layers['counts'] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.layers['lognorm'] = adata.X.copy()

if a.labels == 'final':
    if not (a.label_table and a.label_table_md5):
        sys.exit('--labels final needs --label-table and --label-table-md5')
    got = md5(a.label_table)
    if got != a.label_table_md5:
        sys.exit(f'label table md5 {got} != {a.label_table_md5}')
    t = pd.read_csv(a.label_table, dtype=str, keep_default_na=False)
    assert t['cell_id'].is_unique
    t = t.set_index('cell_id')
    miss = adata.obs_names.difference(t.index)
    if len(miss):
        sys.exit(f'{len(miss)} run cells are not in the label table, e.g. {list(miss[:3])}')
    tt = t.loc[adata.obs_names]
    if not (tt['luminal_object'] == 'yes').all():
        sys.exit('a run cell has luminal_object != yes in the label table')
    # training labels in the table = the run's label column (fact check, not used for scoring)
    if not (tt['scanvi_label'].values == adata.obs[OLD_LABEL_KEY].astype(str).values).all():
        sys.exit('scanvi_label of the table differs from the run\'s label column')
    fin = tt[a.final_column].astype(str)
    if (fin == '').any():
        sys.exit(f'{int((fin == "").sum())} run cells have an empty {a.final_column}')
    adata.obs[LABEL_KEY] = fin.values
    REC.update(label_table=os.path.abspath(a.label_table), label_table_md5=got, final_column=a.final_column,
               final_label_counts_run={k: int(v) for k, v in fin.value_counts().sort_index().items()},
               old_vs_final_crosstab=pd.crosstab(adata.obs[OLD_LABEL_KEY].astype(str).values, fin.values).to_dict())
for key in (BATCH_KEY, LABEL_KEY):
    adata.obs[key] = adata.obs[key].astype(str).fillna('Unknown').astype('category')
log('cells', adata.n_obs, 'genes', adata.n_vars, 'batches', adata.obs[BATCH_KEY].nunique(), 'labels',
    adata.obs[LABEL_KEY].nunique())

# R embeddings as 03
with h5py.File(S0['r_embeddings'], 'r') as f:
    cells = np.array([c.decode() if isinstance(c, bytes) else c for c in f['cells'][:]])
    order = pd.Index(cells).get_indexer(adata.obs_names)
    assert (order >= 0).all(), 'R embeddings do not cover all cells'
    for k in ['X_pca_unintegrated', 'X_umap_unintegrated', 'X_cca', 'X_umap_cca', 'X_rpca', 'X_umap_rpca',
              'X_harmony', 'X_umap_harmony']:
        m = np.asarray(f[k][:])
        if m.shape[0] != len(cells):
            m = m.T
        adata.obsm[k] = m[order]
REC['r_embeddings'] = S0['r_embeddings']
# hvg_var.csv as 03 on a resume
HVG_COLS = ['highly_variable', 'highly_variable_rank', 'means', 'variances', 'variances_norm', 'highly_variable_nbatches']
hv = pd.read_csv(os.path.join(EMB, 'hvg_var.csv'), index_col=0)
assert (hv.index == adata.var_names).all()
for c in HVG_COLS:
    adata.var[c] = hv[c].values
adata.var['highly_variable'] = adata.var['highly_variable'].astype(bool)
# stored Python embeddings as 03 on a resume
for k in ('X_pca_python', 'X_scanorama', 'X_umap_bbknn', 'X_scvi', 'X_scanvi'):
    adata.obsm[k] = np.load(os.path.join(EMB, k + '.npy'))
    assert adata.obsm[k].shape[0] == adata.n_obs
REC['embedding_md5'] = {k: md5(os.path.join(EMB, k + '.npy')) for k in ('X_scanorama', 'X_umap_bbknn', 'X_scvi', 'X_scanvi')}
REC['embedding_md5']['r_embeddings.h5'] = md5(S0['r_embeddings'])
REC['input_md5'] = md5(INPUT)

# ── step 9 as 03 ──
ALL_KEYS = ['X_pca_unintegrated', 'X_cca', 'X_rpca', 'X_harmony', 'X_scanorama', 'X_umap_bbknn', 'X_scvi', 'X_scanvi']
EMBED_KEYS = ALL_KEYS if not a.embeddings else [k for k in ALL_KEYS if k in a.embeddings.split(',')]
assert EMBED_KEYS and (not a.embeddings or len(EMBED_KEYS) == len(a.embeddings.split(','))), a.embeddings
LABEL_MAP = {'X_pca_unintegrated': 'Unintegrated', 'X_cca': 'CCA', 'X_rpca': 'RPCA', 'X_harmony': 'Harmony',
             'X_scanorama': 'Scanorama', 'X_umap_bbknn': 'BBKNN†', 'X_scvi': 'scVI', 'X_scanvi': 'scANVI'}
def rename(df):
    df = df.copy(); df.index = [LABEL_MAP.get(i, i) for i in df.index]; return df
UNLAB = ['Unknown'] if a.labels == 'old' else [v for v in a.unlabelled_values.split(',') if v]
seed_all()
_keep = (~adata.obs[LABEL_KEY].astype(str).isin(UNLAB)).values
bdata = adata[_keep].copy()
for _k in (LABEL_KEY, BATCH_KEY):
    bdata.obs[_k] = bdata.obs[_k].cat.remove_unused_categories()
if set(bdata.obs[LABEL_KEY].astype(str)) & set(UNLAB):
    sys.exit('an unlabeled value is left in the scored set')
REC['scib_scored'] = dict(unlabelled_values=UNLAB, n_cells=int(bdata.n_obs), n_cells_run=int(adata.n_obs),
                          n_unlabelled_run=int((~_keep).sum()),
                          labels={str(k): int(v) for k, v in bdata.obs[LABEL_KEY].value_counts().sort_index().items()},
                          batches=int(bdata.obs[BATCH_KEY].nunique()))
write_csv(pd.DataFrame({'label': bdata.obs[LABEL_KEY].astype(str).values, 'batch': bdata.obs[BATCH_KEY].astype(str).values},
                       index=bdata.obs_names), os.path.join(OUT, 'scib_cells.csv'))
REC['scib_cells_md5'] = md5(os.path.join(OUT, 'scib_cells.csv'))
REC['embedding_keys'] = EMBED_KEYS
log('benchmark cells', json.dumps(REC['scib_scored']))
write_json(REC, os.path.join(OUT, 'rescore_settings.json'))
assert bdata.X is not bdata.layers['counts']
bm = Benchmarker(bdata, batch_key=BATCH_KEY, label_key=LABEL_KEY, embedding_obsm_keys=EMBED_KEYS, n_jobs=1)
bm.benchmark()
write_csv(rename(bm.get_results(min_max_scale=True)), os.path.join(OUT, 'scib_results_minmax_scibmetrics_prefix.csv'))
write_csv(bm._results, os.path.join(OUT, 'scib_results_internal.csv'))
write_csv(rename(bm.get_results(min_max_scale=False)), os.path.join(OUT, 'scib_results_prefix.csv'))
tmp = os.path.join(OUT, 'X_pca_benchmark.tmp.npy'); np.save(tmp, np.asarray(bdata.obsm['X_pca']))
fsync_replace(tmp, os.path.join(OUT, 'X_pca_benchmark.npy'))
log('benchmark\n' + rename(bm.get_results(min_max_scale=False)).to_string())
REC['finished'] = time.strftime('%Y-%m-%d %H:%M:%S')
REC['minutes'] = round((time.time() - T0) / 60, 2)
REC['peak_rss_gb'] = round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 ** 2, 2)
write_json(REC, os.path.join(OUT, 'rescore_settings.json'))
with open(os.path.join(OUT, 'DONE.tmp'), 'w') as fh:
    fh.write(REC['finished'] + '\n')
fsync_replace(os.path.join(OUT, 'DONE.tmp'), os.path.join(OUT, 'DONE'))
log('DONE')
