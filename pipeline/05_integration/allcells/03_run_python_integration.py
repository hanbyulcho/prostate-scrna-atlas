# Code release header: runs six integrations (unintegrated PCA, Harmony, Scanorama, BBKNN, scVI, scANVI), UMAPs and the scib-metrics Benchmarker, for one seed.
# Input: a prepared AnnData object of raw counts (from a preparation script not included). Batch key: Paper_ID_short (the study).
# Output: embeddings/*.npy, tables/scib_results_*.csv, models/ and run_settings.json.
# Defaults for all cells: 2,000 HVGs (seurat_v3, by batch) and 50 PCs. For luminal cells: 3,000 HVGs and 30 PCs.
# In the runs, --label-key was cell_type (the atlas cell types), not the default final_hit1_group_detailed.
# scVI: 30 latent dimensions, 2 layers, up to 400 epochs with early stopping. scANVI: up to 20 epochs. UMAPs use 30 neighbors. Methods: Data integration and batch correction, Benchmarking of integration methods.
"""03_run_python_integration.py - all-cells or luminal Python integration + scib-metrics benchmark.

UMAP sigma fix: umap_sigma_fix.py (in the parent folder) is loaded right after the imports, its guard stops the process unless it passes, and the record goes into
run_settings.json; a resume refuses a folder written without it. No argument, seed or other step changes. BBKNN's graph and layout change, and so do the PCA and Harmony panels and kBET
slightly; PCA, Harmony, Scanorama, scVI, scANVI and the HVG list do not (checked by a test script, not included).

Inputs, seeds and outputs:
  * input: prepared object from 01_prepare_inputs.py (not included; X = raw counts) -> lognorm computed here
    (normalize_total 1e4 + log1p).
  * seeds fixed and recorded (--seed, required), also for scVI and scANVI (scvi.settings.seed).
    --seed 0 gives scanpy PCA/UMAP/neighbors, harmonypy and Scanorama their library default seed 0.
  * scVI and scANVI models are SAVED (model.save, with the HVG list and the setup registry).
  * BBKNN/UMAP on light AnnData objects (a check, not included, gave the same results).
  * per-step .npy outputs (resume); the final object is written only with --write-object.
  * the Benchmarker output is saved before the BBKNN fix; the fix/min-max/figures are step 04.

Thread checks:
  * the thread count changes results (scVI epochs, latents, scib tables), so it is checked and recorded: the run stops
    unless the cores visible to the process (taskset) and OMP/OpenBLAS/MKL/NUMBA thread variables all equal --threads;
    run_settings.json 'thread_settings' records them; the resume guard also refuses a folder written with another
    --threads, scVI/scANVI epoch limit or HVG count.
Run records:
  * run_settings.json records 'run_id' (RUN_ID from the run wrapper, not included) and, in 'thread_settings', 'rule_decision'
    (THREADS_RULE_DECISION: the recorded thread-count rule of the run).
  * the resume guard also compares batch_key, n_pcs, benchmark_n_jobs and the scVI architecture, and a luminal run
    refuses R embeddings whose r_seeds.json thread record differs from --threads.
R-embedding thread check:
  * the R-embedding thread check fails closed: with --r-embeddings, r_seeds.json must exist next to it and record
    threads_env.OMP_NUM_THREADS equal to --threads.
  * thread_settings also records rule_decision_row_md5 (THREADS_RULE_DECISION_MD5 from the run wrapper: md5 of the recorded rule
    that sets the luminal thread count), re-checked by a release step (not included).
All recorded thread variables:
  * the R-embedding thread check compares EVERY thread variable that r_seeds.json records, not only OMP_NUM_THREADS:
    each *_NUM_THREADS / *_MAX_THREADS / *MAXIMUM_THREADS key (02 records OMP, OPENBLAS and RCPP_PARALLEL; MKL, NUMBA
    and any other such key are compared when present) must equal --threads, OMP/OPENBLAS/RCPP_PARALLEL must be present;
    DRYRUN_THREADS, THREADS_PER_SEED_RULE and any other key holding THREAD must equal this process's value of the same
    variable; visible_cores must equal --threads. A record that is not JSON, not an object, or whose threads_env is not
    an object is refused with a message.
Reproducibility and safe writes:
  * all cells: the harmonypy step runs with the BLAS/OpenMP pools limited to 1 thread (threadpoolctl), the only step
    found not bit-reproducible with the same seed at 4 threads; every other step keeps --threads.
  * outputs written to a temp file / folder, fsync'd and renamed (settings, benchmark tables, models, object, DONE).
Seeds and resume:
  * --seed is REQUIRED (no default) and recorded; it reaches every seeded call listed in run_settings.json
    'seed_scheme'. Not seedable here: the BBKNN graph (bbknn 1.6.0 'annoy' search has no seed argument; annoy's own
    fixed default seed applies) and the scib-metrics Benchmarker internals (package defaults).
  * one output folder per seed: the per-step resume (have()) refuses an outdir whose run_settings.json was written
    with another seed or input, so a second seed can never pick up the first seed's .npy files.

Usage: python 03_run_python_integration.py --mode allcells|luminal --input <prepared.h5ad> --outdir <dir>
         [--r-embeddings <r_embeddings.h5>] (luminal) --seed N --threads N [--write-object <out.h5ad>]
         [--scvi-max-epochs 400 --scanvi-max-epochs 20]   (smaller values ONLY for the smoke test)
"""
import argparse, json, os, platform, random, re, resource, sys, time

p = argparse.ArgumentParser()
p.add_argument('--mode', choices=['allcells', 'luminal'], required=True)
p.add_argument('--input', required=True)
p.add_argument('--outdir', required=True)
p.add_argument('--r-embeddings', default=None)
p.add_argument('--batch-key', default='Paper_ID_short')
p.add_argument('--label-key', default=None, help='default: final_hit1_group_detailed (allcells) / etsgroup (luminal)')
p.add_argument('--seed', type=int, required=True, help='required')
p.add_argument('--threads', type=int, required=True, help='required; must equal the visible cores and thread variables')
p.add_argument('--n-hvg', type=int, default=None)
p.add_argument('--n-pcs', type=int, default=None)
p.add_argument('--scvi-max-epochs', type=int, default=400)
p.add_argument('--scanvi-max-epochs', type=int, default=20)
p.add_argument('--scvi-arch', choices=['published', 'scarches'], default='published',
               help="'scarches' = scvi-tools' recommended architecture for reference mapping (NOT the published one)")
p.add_argument('--benchmark-n-jobs', type=int, default=1)
p.add_argument('--skip-benchmark', action='store_true')
p.add_argument('--figure-pdf', action='store_true', help='also write dim_all_methods.pdf (148 MB for all cells)')
p.add_argument('--write-object', default=None)
a = p.parse_args()

# the thread count is a result-determining setting: visible cores and thread variables must agree
_aff = len(os.sched_getaffinity(0))
_tv = {k: os.environ.get(k) for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMBA_NUM_THREADS')}
_bad = {k: v for k, v in _tv.items() if v != str(a.threads)}
if _aff != a.threads or _bad:
    sys.exit(f'--threads {a.threads} but {_aff} visible cores and thread variables {_tv}; run through 08_run_full.sh')
MODE = a.mode
DEF = {  # default settings of each mode
    'allcells': dict(label_key='final_hit1_group_detailed', n_hvg=2000, n_pcs=50, pca_key='X_pca_unintegrated'),
    'luminal': dict(label_key='etsgroup', n_hvg=3000, n_pcs=30, pca_key='X_pca_python'),
}[MODE]
LABEL_KEY = a.label_key or DEF['label_key']
BATCH_KEY = a.batch_key
N_HVG = a.n_hvg or DEF['n_hvg']
N_PCS = a.n_pcs or DEF['n_pcs']
PCA_KEY = DEF['pca_key']
SEED = a.seed

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
# >>> sigma fix: umap_sigma_fix.py (in the parent folder) overrides umap.umap_.smooth_knn_dist on import; guard() stops this process unless all is right
sys.dont_write_bytecode = True; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import umap_sigma_fix
del sys.path[0]
SIGMA_FIX = umap_sigma_fix.guard()
# <<< sigma fix

OUT = os.path.abspath(a.outdir)
EMB = os.path.join(OUT, 'embeddings'); TAB = os.path.join(OUT, 'tables'); FIG = os.path.join(OUT, 'figures')
MOD = os.path.join(OUT, 'models')
# never resume in a folder written with another seed, mode or input (have() only looks at this folder)
_prev = os.path.join(OUT, 'run_settings.json')
if os.path.exists(_prev):
    _s = json.load(open(_prev))['settings']
    _now = dict(seed=SEED, mode=MODE, input=os.path.abspath(a.input), label_key=LABEL_KEY,
                threads=a.threads, n_hvg=N_HVG,   # thread count and HVG count too
                batch_key=BATCH_KEY, n_pcs=N_PCS, benchmark_n_jobs=a.benchmark_n_jobs,
                umap_sigma_fix_md5=SIGMA_FIX.get('module_md5'))   # never resume across the fix
    _sc = _s.get('scvi', {}); _sa = _s.get('scanvi', {})
    if _sc.get('arch') != a.scvi_arch:
        sys.exit(f'{OUT} holds a run with another scVI architecture ({_sc.get("arch")}); use a new --outdir')
    if _sc.get('max_epochs') != a.scvi_max_epochs or _sa.get('max_epochs') != a.scanvi_max_epochs:
        sys.exit(f'{OUT} holds a run with other epoch limits (scVI {_sc.get("max_epochs")}, scANVI {_sa.get("max_epochs")}); '
                 f'use a new --outdir')
    _bad = {k: (_s.get(k), v) for k, v in _now.items() if _s.get(k) != v}
    if _bad:
        sys.exit(f'{OUT} holds a run with other settings {_bad} (previous, now); use a new --outdir per seed')
# R embeddings (luminal) must come from an R step run with the same thread count
if a.r_embeddings:
    _rs = os.path.join(os.path.dirname(os.path.abspath(a.r_embeddings)), 'r_seeds.json')
    _again = 'run the R step again in a new folder'
    try:
        _rj = json.load(open(_rs)) if os.path.exists(_rs) else None
    except ValueError as e:
        sys.exit(f'{a.r_embeddings}: {_rs} is not valid JSON ({e}); {_again}')
    _te = _rj.get('threads_env') if isinstance(_rj, dict) else None
    if _rj is None or _te is None or (isinstance(_te, dict) and _te.get('OMP_NUM_THREADS') is None):
        sys.exit(f'{a.r_embeddings}: no thread record (threads_env.OMP_NUM_THREADS) in {_rs}; the R step cannot be shown '
                 f'to have run with --threads {a.threads}; {_again}')
    if not isinstance(_te, dict):
        sys.exit(f'{a.r_embeddings}: threads_env in {_rs} is a {type(_te).__name__}, not an object; {_again}')
    # every recorded thread variable
    _rx = re.compile(r'.*(?:_NUM_THREADS|_MAX_THREADS|MAXIMUM_THREADS)$')
    _bad = {}
    for _k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'RCPP_PARALLEL_NUM_THREADS'):
        if _k not in _te:
            _bad[_k] = ('not recorded', str(a.threads))
    for _k, _v in _te.items():
        if _rx.fullmatch(str(_k)):
            if str(_v) != str(a.threads):
                _bad[_k] = (_v, str(a.threads))
        elif 'THREAD' in str(_k).upper():
            if str(_v if _v is not None else '') != os.environ.get(str(_k), ''):
                _bad[_k] = (_v, os.environ.get(str(_k), ''))
    _vc = _rj.get('visible_cores')
    if str(_vc) != str(a.threads):
        _bad['visible_cores'] = (_vc, str(a.threads))
    if _bad:
        sys.exit(f'{a.r_embeddings}: the R step thread record ({_rs}) differs from this run: '
                 + '; '.join(f'{k} = {v[0]!r} (here {v[1]!r})' for k, v in sorted(_bad.items()))
                 + f' (this run --threads {a.threads}); {_again}')
for d in (EMB, TAB, FIG, MOD):
    os.makedirs(d, exist_ok=True)
T0 = time.time()
TIMES = {}
def log(*x):
    print(time.strftime('%Y-%m-%d %H:%M:%S'), f'[{MODE} +{(time.time() - T0) / 60:.1f} min]', *x, flush=True)
def peak_gb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 ** 2
class step:
    def __init__(self, name): self.name = name
    def __enter__(self): self.t = time.time(); log('START', self.name); return self
    def __exit__(self, *e):
        TIMES[self.name] = round((time.time() - self.t) / 60, 2)
        log('END', self.name, f'{TIMES[self.name]} min', f'peak RSS {peak_gb():.1f} GB')
def seed_all(s=SEED):
    random.seed(s); np.random.seed(s); torch.manual_seed(s)
def have(*names):
    return all(os.path.exists(os.path.join(EMB, n)) for n in names)
def fsync_replace(tmp, final):
    # fsync the temp file, refuse an empty one, then rename it
    with open(tmp, 'rb+') as fh:
        os.fsync(fh.fileno())
    if os.path.getsize(tmp) == 0:
        sys.exit(f'{tmp} is empty after writing')
    os.replace(tmp, final)
def write_json(obj, final):
    tmp = f'{final}.tmp{os.getpid()}'
    with open(tmp, 'w') as fh:
        json.dump(obj, fh, indent=1)
    fsync_replace(tmp, final)
def write_csv(df, final):
    tmp = f'{final}.tmp{os.getpid()}'; df.to_csv(tmp); fsync_replace(tmp, final)
def save_model(model, final_dir):
    import shutil
    tmp = f'{final_dir}.tmp{os.getpid()}'
    model.save(tmp, save_anndata=False, overwrite=True)
    for f in os.listdir(tmp):
        with open(os.path.join(tmp, f), 'rb+') as fh:
            os.fsync(fh.fileno())
    if os.path.exists(final_dir):
        shutil.rmtree(final_dir)
    os.rename(tmp, final_dir)
def save_npy(name, arr):
    tmp = os.path.join(EMB, name + '.tmp.npy'); np.save(tmp, np.asarray(arr)); fsync_replace(tmp, os.path.join(EMB, name))
def load_npy(name):
    return np.load(os.path.join(EMB, name))

import importlib.metadata as im
VERS = {k: im.version(k) for k in ['scanpy', 'anndata', 'scvi-tools', 'scib-metrics', 'harmonypy', 'scanorama', 'bbknn',
                                   'torch', 'lightning', 'jax', 'jaxlib', 'numpy', 'scipy', 'pandas', 'scikit-learn',
                                   'umap-learn', 'pynndescent', 'numba', 'llvmlite', 'scikit-misc', 'h5py']}
VERS.update(python=platform.python_version(), executable=sys.executable)
SETTINGS = dict(mode=MODE, input=os.path.abspath(a.input), batch_key=BATCH_KEY, label_key=LABEL_KEY, n_hvg=N_HVG,
                n_pcs=N_PCS, hvg='seurat_v3, batch-aware, layer=counts', scale='sc.pp.scale (default)',
                pca_key=PCA_KEY, seed=SEED, threads=a.threads, cpu_affinity=sorted(os.sched_getaffinity(0)),
                thread_settings=dict(threads=a.threads, visible_cores=_aff, env=_tv, torch_threads=torch.get_num_threads(),
                                     rcpp_parallel=os.environ.get('RCPP_PARALLEL_NUM_THREADS'),
                                     rule_threads_per_seed=os.environ.get('THREADS_PER_SEED_RULE'),
                                     rule_decision=os.environ.get('THREADS_RULE_DECISION'),
                                     rule_decision_row_md5=os.environ.get('THREADS_RULE_DECISION_MD5'),
                                     dryrun_thread_override=os.environ.get('DRYRUN_THREADS'),
                                     harmony='1 thread (threadpool_limits) for the harmonypy step only'),
                scvi=dict(n_latent=30, n_layers=2, max_epochs=a.scvi_max_epochs, early_stopping=True,
                          arch=a.scvi_arch, input='raw counts of the HVG subset'),
                scanvi=dict(max_epochs=a.scanvi_max_epochs, n_samples_per_label=100, unlabeled_category='Unknown',
                            init='from_scvi_model'),
                umap_n_neighbors=30, benchmark_n_jobs=a.benchmark_n_jobs, r_embeddings=a.r_embeddings,
                umap_sigma_fix=SIGMA_FIX, umap_sigma_fix_md5=SIGMA_FIX.get('module_md5'),
                seed_scheme=dict(  # where --seed goes (recorded in run_settings.json)
                    seed_all='python random, numpy, torch.manual_seed before every stochastic step',
                    pca='sc.pp.pca random_state', harmonypy='run_harmony random_state (all cells)',
                    scanorama='integrate_scanpy seed', bbknn_umap='sc.tl.umap random_state',
                    harmony_threads='1 (threadpool_limits(1) for the harmonypy step only)',
                    scvi_scanvi='scvi.settings.seed', umaps='sc.pp.neighbors / sc.tl.umap random_state',
                    not_seeded_here=['BBKNN graph: bbknn 1.6.0 annoy search has no seed argument (annoy default seed)',
                                     'scib-metrics Benchmarker internals: package defaults'],
                    luminal_R='02_r_seurat_luminal.R --seed (r_seeds.json)'),
                run_id=os.environ.get('RUN_ID'),   # read by a run-registry script (not included)
                started=time.strftime('%Y-%m-%d %H:%M:%S'))
write_json(dict(settings=SETTINGS, versions=VERS), os.path.join(OUT, 'run_settings.json'))
log('settings', json.dumps(SETTINGS)); log('versions', json.dumps(VERS))

# ── Step 1: load counts, lognorm ────────────────────────────────────────────────────────────────────────────────
with step('1_load_lognorm'):
    adata = ad.read_h5ad(a.input)
    assert sp.issparse(adata.X)
    d = adata.X.data
    assert np.all(d == np.floor(d)), 'X must hold integer raw counts (01_prepare_inputs.py output)'
    adata.layers['counts'] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.layers['lognorm'] = adata.X.copy()
    for key in (BATCH_KEY, LABEL_KEY):
        adata.obs[key] = adata.obs[key].astype(str).fillna('Unknown').astype('category')
    log('cells', adata.n_obs, 'genes', adata.n_vars, 'batches', adata.obs[BATCH_KEY].nunique(),
        'labels', adata.obs[LABEL_KEY].nunique(), 'Unknown labels', int((adata.obs[LABEL_KEY] == 'Unknown').sum()))

if MODE == 'luminal':
    assert a.r_embeddings, '--r-embeddings (02_r_seurat_luminal.R output) is required for --mode luminal'
    with h5py.File(a.r_embeddings, 'r') as f:
        cells = np.array([c.decode() if isinstance(c, bytes) else c for c in f['cells'][:]])
        order = pd.Index(cells).get_indexer(adata.obs_names)
        assert (order >= 0).all(), 'R embeddings do not cover all cells'
        for k in ['X_pca_unintegrated', 'X_umap_unintegrated', 'X_cca', 'X_umap_cca', 'X_rpca', 'X_umap_rpca',
                  'X_harmony', 'X_umap_harmony']:
            m = np.asarray(f[k][:])
            if m.shape[0] != len(cells):
                m = m.T
            adata.obsm[k] = m[order]
    log('R embeddings', {k: adata.obsm[k].shape for k in adata.obsm})

# ── Step 2: HVG (seurat_v3 on counts, batch-aware) + scaled PCA ──────────────────────────────────────────────────
HVG_COLS = ['highly_variable', 'highly_variable_rank', 'means', 'variances', 'variances_norm', 'highly_variable_nbatches']
with step('2_hvg_pca'):
    if have('hvg_var.csv'):
        hv = pd.read_csv(os.path.join(EMB, 'hvg_var.csv'), index_col=0)
        assert (hv.index == adata.var_names).all()
        for c in HVG_COLS:
            adata.var[c] = hv[c].values
        adata.var['highly_variable'] = adata.var['highly_variable'].astype(bool)
    else:
        sc.pp.highly_variable_genes(adata, n_top_genes=N_HVG, batch_key=BATCH_KEY, flavor='seurat_v3', layer='counts')
        write_csv(adata.var[HVG_COLS], os.path.join(EMB, 'hvg_var.csv'))
    log('HVG', int(adata.var['highly_variable'].sum()))
    adata_hvg = adata[:, adata.var['highly_variable'].values].copy()
    del adata_hvg.layers['counts']
    sc.pp.scale(adata_hvg)
    if have(PCA_KEY + '.npy'):
        adata.obsm[PCA_KEY] = load_npy(PCA_KEY + '.npy')
    else:
        seed_all()
        sc.pp.pca(adata_hvg, n_comps=N_PCS, random_state=SEED)
        adata.obsm[PCA_KEY] = adata_hvg.obsm['X_pca'].copy()
        save_npy(PCA_KEY + '.npy', adata.obsm[PCA_KEY])
    log(PCA_KEY, adata.obsm[PCA_KEY].shape)

# ── Step 3: Harmony (harmonypy; all cells only - the luminal Harmony is Seurat/R) ────────────────────────────────
if MODE == 'allcells':
    with step('3_harmony'):
        if have('X_harmony.npy'):
            adata.obsm['X_harmony'] = load_npy('X_harmony.npy')
        else:
            seed_all()
            # harmonypy is bit-reproducible only single-threaded, so the BLAS and OpenMP
            # pools (numpy, sklearn KMeans) are limited to 1 thread for this step only; torch threads are left unchanged
            from threadpoolctl import threadpool_limits, threadpool_info
            with threadpool_limits(limits=1):
                HARMONY_POOLS = [dict(api=x.get('internal_api'), threads=x.get('num_threads')) for x in threadpool_info()]
                ho = hm.run_harmony(adata.obsm[PCA_KEY], adata.obs, BATCH_KEY, random_state=SEED)
            log('harmonypy thread pools during the step:', json.dumps(HARMONY_POOLS))
            SETTINGS['harmony_threadpools'] = HARMONY_POOLS
            adata.obsm['X_harmony'] = np.asarray(ho.Z_corr.T)
            save_npy('X_harmony.npy', adata.obsm['X_harmony'])

# ── Step 4: Scanorama (per-batch HVG subsets, batches sorted) ────────────────────────────────────────────────────
with step('4_scanorama'):
    if have('X_scanorama.npy'):
        adata.obsm['X_scanorama'] = load_npy('X_scanorama.npy')
    else:
        seed_all()
        batches_sorted = sorted(adata_hvg.obs[BATCH_KEY].unique())
        adatas_list = [adata_hvg[adata_hvg.obs[BATCH_KEY] == b].copy() for b in batches_sorted]
        scanorama.integrate_scanpy(adatas_list, dimred=N_PCS, seed=SEED)
        emb = pd.concat([pd.DataFrame(x.obsm['X_scanorama'], index=x.obs_names) for x in adatas_list])
        adata.obsm['X_scanorama'] = emb.loc[adata.obs_names].values
        save_npy('X_scanorama.npy', adata.obsm['X_scanorama'])
        del adatas_list, emb
del adata_hvg

# ── Step 5: BBKNN on the PCA + UMAP of the BBKNN graph (scib proxy embedding) ────────────────────────────────────
with step('5_bbknn'):
    if have('X_umap_bbknn.npy', 'bbknn_connectivities.npz'):
        adata.obsm['X_umap_bbknn'] = load_npy('X_umap_bbknn.npy')
        adata.obsp['bbknn_connectivities'] = sp.load_npz(os.path.join(EMB, 'bbknn_connectivities.npz'))
    else:
        seed_all()
        b = ad.AnnData(obs=adata.obs[[BATCH_KEY]].copy())
        b.obsm['X_pca'] = adata.obsm[PCA_KEY].copy()
        bbknn.bbknn(b, batch_key=BATCH_KEY, n_pcs=N_PCS)
        sc.tl.umap(b, random_state=SEED)
        adata.obsm['X_umap_bbknn'] = b.obsm['X_umap'].copy()
        adata.obsp['bbknn_connectivities'] = b.obsp['connectivities'].tocsr()
        sp.save_npz(os.path.join(EMB, 'bbknn_connectivities.npz'), adata.obsp['bbknn_connectivities'])
        save_npy('X_umap_bbknn.npy', adata.obsm['X_umap_bbknn'])
        del b

# ── Step 6: scVI + scANVI on the HVG raw counts; models SAVED ───────────────────────────────────────────────────
SCVI_DIR, SCANVI_DIR = os.path.join(MOD, 'scvi'), os.path.join(MOD, 'scanvi')
ARCH = {} if a.scvi_arch == 'published' else dict(use_layer_norm='both', use_batch_norm='none',
                                                    encode_covariates=True, dropout_rate=0.2)
def make_scvi_adata():
    hv = adata.var['highly_variable'].values
    s = ad.AnnData(X=adata.layers['counts'][:, np.where(hv)[0]].copy(), obs=adata.obs[[BATCH_KEY, LABEL_KEY]].copy(),
                   var=pd.DataFrame(index=adata.var_names[hv]))
    scvi.model.SCVI.setup_anndata(s, batch_key=BATCH_KEY)
    return s
model_scvi = None
with step('6a_scvi'):
    if have('X_scvi.npy') and os.path.exists(os.path.join(SCVI_DIR, 'model.pt')):
        adata.obsm['X_scvi'] = load_npy('X_scvi.npy')
    else:
        s = make_scvi_adata()
        seed_all(); scvi.settings.seed = SEED
        model_scvi = scvi.model.SCVI(s, n_latent=30, n_layers=2, **ARCH)
        model_scvi.train(max_epochs=a.scvi_max_epochs, early_stopping=True, datasplitter_kwargs={'num_workers': 0})
        adata.obsm['X_scvi'] = model_scvi.get_latent_representation()
        save_model(model_scvi, SCVI_DIR)
        write_csv(pd.concat({k: v for k, v in model_scvi.history.items()}, axis=1), os.path.join(MOD, 'scvi_history.csv'))
        save_npy('X_scvi.npy', adata.obsm['X_scvi'])
        log('scVI epochs trained', len(model_scvi.history['elbo_train']))
with step('6b_scanvi'):
    if have('X_scanvi.npy') and os.path.exists(os.path.join(SCANVI_DIR, 'model.pt')):
        adata.obsm['X_scanvi'] = load_npy('X_scanvi.npy')
    else:
        s = make_scvi_adata()
        if model_scvi is None:
            model_scvi = scvi.model.SCVI.load(SCVI_DIR, adata=s)
        s.obs['celltype_scanvi'] = adata.obs[LABEL_KEY].astype(str).values
        seed_all(); scvi.settings.seed = SEED
        model_scanvi = scvi.model.SCANVI.from_scvi_model(model_scvi, adata=s, labels_key='celltype_scanvi',
                                                         unlabeled_category='Unknown')
        model_scanvi.train(max_epochs=a.scanvi_max_epochs, n_samples_per_label=100,
                           datasplitter_kwargs={'num_workers': 0})
        adata.obsm['X_scanvi'] = model_scanvi.get_latent_representation()
        save_model(model_scanvi, SCANVI_DIR)
        write_csv(pd.concat({k: v for k, v in model_scanvi.history.items()}, axis=1), os.path.join(MOD, 'scanvi_history.csv'))
        write_csv(pd.DataFrame({'label': s.obs['celltype_scanvi'].values, 'scanvi_prediction': model_scanvi.predict()},
                               index=s.obs_names), os.path.join(MOD, 'scanvi_predictions.csv'))
        save_npy('X_scanvi.npy', adata.obsm['X_scanvi'])
        del model_scanvi
model_scvi = None
write_json(dict(seed=SEED, scvi_settings_seed=SEED, torch_manual_seed=SEED, numpy_seed=SEED, python_random_seed=SEED,
                scvi_dir=SCVI_DIR, scanvi_dir=SCANVI_DIR, hvg_genes=list(adata.var_names[adata.var['highly_variable']]),
                batch_categories=list(adata.obs[BATCH_KEY].cat.categories),
                label_categories=list(adata.obs[LABEL_KEY].cat.categories)),
           os.path.join(MOD, 'model_manifest.json'))

# ── Step 7: UMAPs (n_neighbors=30) ─────────────────────────────────────────────────────────────────────────────
UMAPS = ([('X_pca_unintegrated', 'X_umap_unintegrated'), ('X_harmony', 'X_umap_harmony')] if MODE == 'allcells' else []) + \
        [('X_scanorama', 'X_umap_scanorama'), ('X_scvi', 'X_umap_scvi'), ('X_scanvi', 'X_umap_scanvi')]
with step('7_umaps'):
    for ek, uk in UMAPS:
        if have(uk + '.npy'):
            adata.obsm[uk] = load_npy(uk + '.npy'); continue
        seed_all()
        t = ad.AnnData(obs=pd.DataFrame(index=adata.obs_names)); t.obsm['X_emb'] = adata.obsm[ek]
        sc.pp.neighbors(t, use_rep='X_emb', n_neighbors=30, random_state=SEED)
        sc.tl.umap(t, random_state=SEED)
        adata.obsm[uk] = t.obsm['X_umap'].copy(); save_npy(uk + '.npy', adata.obsm[uk])
        log(uk, 'done')

# ── Step 8: overview figure (dim_all_methods) ──────────────────────────────────────
with step('8_figure'):
    try:
        paper_col = {"P01": "#009CD2", "P02": "#e6194B", "P03": "#3cb44b", "P05": "#ffe119", "P06": "#4363d8",
                     "P07": "#f58231", "P08": "#911eb4", "P09": "#46f0f0", "P10": "#2b1700", "P11": "#f032e6",
                     "P12": "#bcf60c", "P14": "#fabebe", "P15": "#008080", "P16": "#e6beff", "P17": "#9A6324",
                     "P19": "#FFFF33", "P21": "#999999"}
        label_col = {"Luminal": "#E31A1C", "Club": "#E2ADB2", "Basal": "#1976D2", "Hillock": "#87CEFA", "Ionocyte": "#4A1486",
                     "T cell": "#3B7A57", "B cell": "#4DAF4A", "Plasma": "#008B8B", "NK cell": "#00F0A8",
                     "Granulocytes": "#8FBC8F", "Macrophage": "#00AB66", "Monocyte": "#CEFF00", "Neutrophils": "#997950",
                     "Dendritic": "#5D4037", "Endothelial cell": "#FFD800", "SMC": "#FD8D3C", "Fibroblast": "#F0E68C",
                     "Schwann": "#696969", "N": "#71A45E", "Benign": "#87CEFA", "BPH": "#0047ab", "nonETS": "#F781BF",
                     "ERG": "#984EA3", "ETV": "#FF7F00", "AR": "#E41A1C", "NEPC": "#A65628", "DN": "#4DAF4A"}
        methods = ([('X_umap_unintegrated', 'Unintegrated (PCA)')] +
                   ([('X_umap_harmony', 'Harmony')] if MODE == 'allcells' else
                    [('X_umap_cca', 'CCA (Seurat)'), ('X_umap_rpca', 'RPCA (Seurat)'), ('X_umap_harmony', 'Harmony (Seurat)')]) +
                   [('X_umap_scanorama', 'Scanorama'), ('X_umap_bbknn', 'BBKNN'), ('X_umap_scvi', 'scVI'), ('X_umap_scanvi', 'scANVI')])
        q = ad.AnnData(obs=adata.obs[[BATCH_KEY, LABEL_KEY]].copy())
        q.obs['Paper_ID_short2'] = q.obs[BATCH_KEY].astype(str).str.extract(r'(P\d+)')[0].astype('category')
        for k, _ in methods:
            q.obsm[k] = adata.obsm[k]
        fig, axes = plt.subplots(len(methods), 2, figsize=(14, 5 * len(methods)))
        for r, (k, title) in enumerate(methods):
            sc.pl.embedding(q, basis=k, color='Paper_ID_short2', ax=axes[r, 0], show=False, frameon=False, s=3,
                            palette=[paper_col.get(c, '#AAAAAA') for c in q.obs['Paper_ID_short2'].cat.categories],
                            title=f'{title}\nBatch (Paper)')
            sc.pl.embedding(q, basis=k, color=LABEL_KEY, ax=axes[r, 1], show=False, frameon=False, s=3,
                            palette=[label_col.get(c, '#AAAAAA') for c in q.obs[LABEL_KEY].cat.categories],
                            title=f'{title}\n{LABEL_KEY}')
        plt.tight_layout()
        fig.savefig(os.path.join(FIG, 'dim_all_methods.png'), dpi=100, bbox_inches='tight')
        if a.figure_pdf:
            fig.savefig(os.path.join(FIG, 'dim_all_methods.pdf'), bbox_inches='tight')
        plt.close('all')
    except Exception as e:
        log('figure failed (non-fatal):', repr(e))

# ── Step 9: scib-metrics Benchmarker (defaults; pre-integration PCA computed by the Benchmarker on lognorm X) ────
EMBED_KEYS = (['X_pca_unintegrated', 'X_harmony', 'X_scanorama', 'X_umap_bbknn', 'X_scvi', 'X_scanvi'] if MODE == 'allcells'
              else ['X_pca_unintegrated', 'X_cca', 'X_rpca', 'X_harmony', 'X_scanorama', 'X_umap_bbknn', 'X_scvi', 'X_scanvi'])
LABEL_MAP = {'X_pca_unintegrated': 'Unintegrated', 'X_cca': 'CCA', 'X_rpca': 'RPCA', 'X_harmony': 'Harmony',
             'X_scanorama': 'Scanorama', 'X_umap_bbknn': 'BBKNN†', 'X_scvi': 'scVI', 'X_scanvi': 'scANVI'}
def rename(df):
    df = df.copy(); df.index = [LABEL_MAP.get(i, i) for i in df.index]; return df
if not a.skip_benchmark:
    with step('9_benchmark'):
        if os.path.exists(os.path.join(TAB, 'scib_results_prefix.csv')):
            log('benchmark exists, skipped')
        else:
            seed_all()
            assert adata.X is not adata.layers['counts']
            bm = Benchmarker(adata, batch_key=BATCH_KEY, label_key=LABEL_KEY, embedding_obsm_keys=EMBED_KEYS,
                             n_jobs=a.benchmark_n_jobs)
            bm.benchmark()
            write_csv(rename(bm.get_results(min_max_scale=True)), os.path.join(TAB, 'scib_results_minmax_scibmetrics_prefix.csv'))
            write_csv(bm._results, os.path.join(TAB, 'scib_results_internal.csv'))
            write_csv(rename(bm.get_results(min_max_scale=False)), os.path.join(TAB, 'scib_results_prefix.csv'))  # last: resume marker
            save_npy('X_pca_benchmark.npy', adata.obsm['X_pca'])
            log('benchmark\n' + rename(bm.get_results(min_max_scale=False)).to_string())

# ── Step 10: final object (optional) ──────────────────────────────────────────────────────────────────────────
if a.write_object:
    with step('10_write_object'):
        adata.X = adata.layers['lognorm']
        adata.uns['reintegration'] = json.dumps(dict(settings=SETTINGS, versions=VERS, step_minutes=TIMES))
        _tmp = a.write_object[:-5] + f'.tmp{os.getpid()}.h5ad' if a.write_object.endswith('.h5ad') else a.write_object + f'.tmp{os.getpid()}'
        adata.write_h5ad(_tmp); fsync_replace(_tmp, a.write_object)
        log('wrote', a.write_object)

SETTINGS['finished'] = time.strftime('%Y-%m-%d %H:%M:%S')
write_json(dict(settings=SETTINGS, versions=VERS, step_minutes=TIMES, peak_rss_gb=round(peak_gb(), 2)),
           os.path.join(OUT, 'run_settings.json'))
with open(os.path.join(OUT, 'DONE.tmp'), 'w') as fh:
    fh.write(SETTINGS['finished'] + '\n')
fsync_replace(os.path.join(OUT, 'DONE.tmp'), os.path.join(OUT, 'DONE'))
log('ALL DONE', json.dumps(TIMES))
