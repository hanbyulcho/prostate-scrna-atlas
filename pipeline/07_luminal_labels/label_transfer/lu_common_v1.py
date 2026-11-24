# Code release header: shared functions of S1_soft_knn_primary_v1.py (label transfer for luminal cells).
# soft_predict() reloads a saved scANVI model and returns its soft class probabilities, with no training.
# knn_transfer() runs the scArches 0.6.1 weighted kNN (k = 50, Euclidean, brute force) from the installed package.
# Uncertainty is 1 minus the kernel-weighted share of the winning label among the k neighbors. Methods: Data integration and batch correction, Luminal subgroup annotation.
#
# The vectorized check in knn_transfer() recomputes the kernel weights in the same way as scarches/utils/knn.py of scArches 0.6.1.
# Its license notice, from the installed scArches 0.6.1 distribution, follows.
#
# BSD 3-Clause License
#
# Copyright (c) 2022, Mohammad Lotfollahi, Theislab
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# lu_common_v1.py -- shared functions of the luminal label-uncertainty check.
# Nothing here trains a model or changes a label.
#
# soft_predict(): reload a SAVED scANVI model on the stored run input and call scvi-tools
#   SCANVI.predict(soft=True) (scvi-tools 1.4.0.post1, scvi/model/base/_training_mixin.py l.196: softmax over the
#   classifier logits; use_posterior_mean=True, the default, as in the run's own model_scanvi.predict()).
#   The AnnData is built exactly as the run's saved-model check (07_check_saved_models.py, not included) builds it.
#   Checks: reloaded latent equals the stored embeddings/X_scanvi.npy (max |diff| <= 1e-4, the tolerance of that check) and the
#   argmax of the soft output equals the stored hard prediction (python/models/scanvi_predictions.csv) for every cell.
# knn_transfer(): scArches 0.6.1 scarches/utils/knn.py weighted_knn_trainer + weighted_knn_transfer (L. Sikkema's
#   code, used for the HLCA), called unchanged: k = 50 (the function default), Euclidean, brute force; reference = the
#   run's labeled cells (training label not 'Unknown') in the stored X_scanvi of the same seed; query = the cells given.
#   uncertainty = 1 - (kernel-weighted share of the winning label among the 50 neighbors) ("package" mode).
import os, json, hashlib
import numpy as np, pandas as pd

CLASSES = ['AR', 'BPH', 'Benign', 'ERG', 'ETV', 'NEPC', 'nonETS']
SITE = '<PY311_ENV>/lib/python3.11/site-packages'
KNN_PY = f'{SITE}/scarches/utils/knn.py'


def md5(p, bs=1 << 22):
    h = hashlib.md5()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(bs), b''):
            h.update(b)
    return h.hexdigest()


def write_csv_safe(df, path, header_line=None, compress=False):
    """temp file -> flush+fsync -> read back and compare shape/columns -> os.replace."""
    tmp = path + f'.tmp{os.getpid()}'
    if compress:
        import gzip
        with gzip.open(tmp, 'wt', newline='') as fh:
            if header_line:
                fh.write('"' + header_line.replace('\n', ' ').replace('"', '""') + '"\n')
            df.to_csv(fh, index=False)
        with open(tmp, 'rb+') as fh:
            os.fsync(fh.fileno())
    else:
        with open(tmp, 'w', newline='') as fh:
            if header_line:
                fh.write('"' + header_line.replace('\n', ' ').replace('"', '""') + '"\n')
            df.to_csv(fh, index=False)
            fh.flush(); os.fsync(fh.fileno())
    back = pd.read_csv(tmp, skiprows=1 if header_line else 0, low_memory=False, keep_default_na=False, dtype=str,
                       compression='gzip' if compress else None)
    assert back.shape == df.shape, (path, back.shape, df.shape)
    assert list(back.columns) == [str(c) for c in df.columns], path
    os.replace(tmp, path)
    print('wrote', path, df.shape, flush=True)


def write_json_safe(obj, path):
    tmp = path + f'.tmp{os.getpid()}'
    with open(tmp, 'w') as fh:
        json.dump(obj, fh, indent=1, default=str); fh.flush(); os.fsync(fh.fileno())
    json.load(open(tmp))
    os.replace(tmp, path)


def soft_predict(seed_dir, tol=1e-4, threads=8):
    """seed_dir = <run>/luminal_seed<s>. Returns (probs DataFrame indexed by cell, check dict)."""
    import anndata as ad
    import torch, scvi
    torch.set_num_threads(threads)
    run = os.path.join(seed_dir, 'python')
    man = json.load(open(os.path.join(run, 'models', 'model_manifest.json')))
    st = json.load(open(os.path.join(run, 'run_settings.json')))['settings']
    B, Lk = st['batch_key'], st['label_key']
    inp = st['input']
    assert os.path.realpath(inp) == os.path.realpath(os.path.join(seed_dir, 'input.h5ad')), inp
    a_in = ad.read_h5ad(inp)
    q = a_in[:, man['hvg_genes']].copy()
    del a_in
    q.obs = q.obs[[B, Lk]].copy()
    q.obs['celltype_scanvi'] = q.obs[Lk].astype(str).values
    m = scvi.model.SCANVI.load(man['scanvi_dir'], adata=q)
    z = m.get_latent_representation()
    ref = np.load(os.path.join(run, 'embeddings', 'X_scanvi.npy'))
    dz = float(np.abs(z - ref).max())
    assert dz <= tol, f'reloaded scANVI latent differs from the stored X_scanvi: {dz}'
    P = m.predict(soft=True)
    assert isinstance(P, pd.DataFrame), type(P)
    P.index = q.obs_names
    cols = [str(c) for c in P.columns]
    assert sorted(cols) == sorted(CLASSES), cols
    P = P[CLASSES].astype('float64')
    rs = P.sum(axis=1).values
    stored = pd.read_csv(os.path.join(run, 'models', 'scanvi_predictions.csv'), index_col=0)
    assert list(stored.index) == list(P.index), 'cell order differs from scanvi_predictions.csv'
    hard = P.idxmax(axis=1).values
    n_eq = int((hard == stored['scanvi_prediction'].astype(str).values).sum())
    chk = dict(seed_dir=seed_dir, input=inp, scanvi_dir=man['scanvi_dir'], scanvi_model_md5=md5(os.path.join(man['scanvi_dir'], 'model.pt')),
               n_cells=int(P.shape[0]), latent_max_abs_diff_vs_stored=dz, row_sum_min=float(rs.min()), row_sum_max=float(rs.max()),
               argmax_equals_stored_hard_prediction=n_eq, label_equals_stored_label=int((q.obs['celltype_scanvi'].values == stored['label'].astype(str).values).sum()),
               scvi_tools=scvi.__version__, use_posterior_mean=True)
    assert n_eq == P.shape[0], f'argmax of soft output differs from stored hard prediction in {P.shape[0]-n_eq} cells'
    out = P.copy()
    out.insert(0, 'training_label', q.obs['celltype_scanvi'].values)
    out.insert(0, 'cell_id', P.index.values)
    out['p_max'] = P.max(axis=1).values
    out['argmax'] = hard
    return out.reset_index(drop=True), chk


def knn_transfer(emb, labels, query_mask, cell_ids, k=50):
    """scArches 0.6.1 weighted kNN, unchanged. emb: (n, d) array; labels: training labels with 'Unknown' for unlabeled."""
    import anndata as ad, importlib.util
    # the scarches package __init__ fails under anndata 0.12.4 ('cannot import name read'), so the module file
    # scarches/utils/knn.py of the installed scArches 0.6.1 is loaded on its own; its code is used unchanged.
    spec = importlib.util.spec_from_file_location('scarches_knn_0_6_1', KNN_PY)
    kmod = importlib.util.module_from_spec(spec); spec.loader.exec_module(kmod)
    weighted_knn_trainer, weighted_knn_transfer = kmod.weighted_knn_trainer, kmod.weighted_knn_transfer
    labels = np.asarray(labels).astype(str)
    ref_mask = labels != 'Unknown'
    ref = ad.AnnData(X=np.zeros((int(ref_mask.sum()), 1), dtype='float32'))
    ref.obsm['emb'] = np.asarray(emb[ref_mask], dtype='float64')
    ref.obs_names = [str(x) for x in np.asarray(cell_ids)[ref_mask]]
    qry = ad.AnnData(X=np.zeros((int(query_mask.sum()), 1), dtype='float32'))
    qry.obsm['emb'] = np.asarray(emb[query_mask], dtype='float64')
    qry.obs_names = [str(x) for x in np.asarray(cell_ids)[query_mask]]
    ref_obs = pd.DataFrame({'trainlabel': labels[ref_mask]}, index=ref.obs_names)
    knn = weighted_knn_trainer(ref, 'emb', n_neighbors=k)
    pred, unc = weighted_knn_transfer(qry, 'emb', ref_obs, 'trainlabel', knn, threshold=1, pred_unknown=False, mode='package')
    out = pd.DataFrame({'cell_id': qry.obs_names, 'knn_label': pred['trainlabel'].astype(str).values,
                        'knn_uncertainty': unc['trainlabel'].astype(float).values})
    assert out['knn_label'].isin(CLASSES).all() and out['knn_uncertainty'].notna().all()
    # independent recomputation of the same formula (vectorized), to prove the chained pandas assignment inside the
    # package wrote every row (pandas 2.3.3, no copy-on-write)
    d, idx = knn.kneighbors(X=qry.obsm['emb'])
    s = np.std(d, axis=1); s = ((2.0 / s) ** 2).reshape(-1, 1)
    w = np.exp(-np.true_divide(d, s)); w = w / w.sum(axis=1, keepdims=True)
    lab = labels[ref_mask][idx]
    best = np.zeros(len(w)); bestlab = np.empty(len(w), dtype=object)
    for c in sorted(set(lab.ravel())):
        sc_ = (w * (lab == c)).sum(axis=1)
        upd = sc_ > best
        best[upd] = sc_[upd]; bestlab[upd] = c
    rec_unc = np.maximum(1 - best, 0)
    chk = dict(scarches='0.6.1 (dist-info)', knn_py=KNN_PY, knn_py_md5=md5(KNN_PY), k=k, n_reference=int(ref_mask.sum()), n_query=int(query_mask.sum()),
               max_abs_diff_uncertainty_vs_recomputed=float(np.abs(rec_unc - out['knn_uncertainty'].values).max()),
               label_equal_recomputed=int((bestlab == out['knn_label'].values).sum()))
    return out, chk
