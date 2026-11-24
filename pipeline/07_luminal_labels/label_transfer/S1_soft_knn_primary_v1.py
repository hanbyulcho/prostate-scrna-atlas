# Code release header: scArches weighted kNN label transfer (k = 50) in luminal scANVI space (X_scanvi), seeds 0 and 1, with soft scANVI probabilities.
# Input: each seed's stored run (scANVI model, input object, X_scanvi, predictions). Reference: cells with a training class. Query: cells without one ('Unknown', released as 'unlabeled').
# Output: primary_seed<s>_soft_v1.csv.gz and primary_seed<s>_knn_v1.csv.gz for each seed s.
# The output was curator evidence for four units. Cells with uncertainty above 0.2 (a cutoff chosen post hoc) counted as unknown in each unit's majority.
# A unit label needed both seeds to agree and the same label on a second luminal integration that left the unit out (script not included). Methods: Data integration and batch correction, Luminal subgroup annotation.
# S1_soft_knn_primary_v1.py -- step 1 of the luminal label-uncertainty check.
# For the primary luminal integration run, seeds 0 and 1:
#   (a) soft scANVI probabilities from the SAVED model (lu_common_v1.soft_predict; no training);
#   (b) scArches weighted-kNN label transfer with uncertainty on the stored X_scanvi (lu_common_v1.knn_transfer;
#       reference = labeled cells, query = every 'Unknown' cell of the run).
# Writes to the WORK folder only.
import os, sys, json, time
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lu_common_v1 import soft_predict, knn_transfer, write_csv_safe, write_json_safe, md5

REV = '<PROJECT_DIR>'
RUN = f'{REV}/reintegration_v18/runs/20261107_luminal_primary'
WORK = '<SCRATCH_DIR>/lum_uncertainty_v1/work'
PIN = {  # stored outputs of the run, with their md5
    'pred_s0': (f'{RUN}/luminal_seed0/python/models/scanvi_predictions.csv', '0950ef005c3c6880ceacc86b41cd614e'),
    'pred_s1': (f'{RUN}/luminal_seed1/python/models/scanvi_predictions.csv', 'c4b64cfe767a9110fd8fcd29b035a9c9'),
}
TH = int(os.environ.get('LU_THREADS', '8'))
log = dict(started=time.strftime('%Y-%m-%d %H:%M:%S %Z'), run=RUN, threads=TH, seeds={})
for k, (p, m) in PIN.items():
    got = md5(p); assert got == m, (k, got, m)
log['pins'] = {k: v[1] for k, v in PIN.items()}
for s in (0, 1):
    sd = f'{RUN}/luminal_seed{s}'
    t0 = time.time()
    P, chk = soft_predict(sd, threads=TH)
    chk['minutes_soft'] = round((time.time() - t0) / 60, 2)
    write_csv_safe(P, f'{WORK}/primary_seed{s}_soft_v1.csv.gz', compress=True)
    t1 = time.time()
    emb = np.load(f'{sd}/python/embeddings/X_scanvi.npy')
    assert emb.shape[0] == len(P)
    qmask = (P['training_label'].values == 'Unknown')
    K, kchk = knn_transfer(emb, P['training_label'].values, qmask, P['cell_id'].values, k=50)
    kchk['minutes_knn'] = round((time.time() - t1) / 60, 2)
    kchk['X_scanvi_md5'] = md5(f'{sd}/python/embeddings/X_scanvi.npy')
    write_csv_safe(K, f'{WORK}/primary_seed{s}_knn_v1.csv.gz', compress=True)
    log['seeds'][s] = dict(soft=chk, knn=kchk)
    print(json.dumps(log['seeds'][s], indent=1), flush=True)
log['finished'] = time.strftime('%Y-%m-%d %H:%M:%S %Z')
write_json_safe(log, f'{WORK}/S1_log_v1.json')
print('S1 done', flush=True)
