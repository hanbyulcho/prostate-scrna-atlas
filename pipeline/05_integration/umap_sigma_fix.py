# Code release header: umap-learn 0.5.9.post2 with numba 0.62.1 gives sigma = inf in UMAP's fuzzy simplicial set for cells whose correct sigma is above 1.
# Importing this module recompiles umap-learn's smooth_knn_dist without fastmath, so the UMAPs, BBKNN graphs and kBET scores stay correct.
# Both 03_run_python_integration.py scripts and R1_rescore_scib.py load it, and its guard() stops the run unless the patch is in place and correct on test inputs.
"""umap_sigma_fix.py - replaces umap.umap_.smooth_knn_dist with the same function compiled without fastmath.

The fault: in the human_sc venv (numba 0.62.1, llvmlite 0.45.1, umap-learn 0.5.9.post2) the compiled umap.umap_.smooth_knn_dist returns
sigma = inf for any cell whose correct sigma is above 1, so all that cell's neighbors get weight exactly 1. The function is compiled with
fastmath=True, and its search doubles mid while hi == NPY_INFINITY. A likely cause (not verified): fastmath lets the
compiler assume there are no infinities.

The fix: importing this module replaces umap.umap_.smooth_knn_dist with the SAME source and the same locals, compiled with fastmath=False.
No argument changes. Nothing is pinned to other package versions.

Why this reaches every caller: fuzzy_simplicial_set is a plain Python function in umap.umap_ that looks up smooth_knn_dist in its module's
globals at call time. So names bound at import (from umap.umap_ import fuzzy_simplicial_set, in bbknn.matrix and in
scib_metrics.nearest_neighbors._dataclass; at call time in scanpy.neighbors._connectivity) still reach the override, whichever is imported first.
The override lives in one process only: a child process that reaches fuzzy_simplicial_set must import this module itself.
03_run_python_integration.py starts no child process that does (checked by a test script, not included).

Use, right after the imports of a process, before any first use:
    sys.dont_write_bytecode = True; sys.path.insert(0, <this folder>)     # python -I leaves the script's folder out of sys.path
    import umap_sigma_fix
    del sys.path[0]
    SIGMA_FIX = umap_sigma_fix.guard()      # stops the process unless all is right; goes into the run settings
There is no switch to turn this off."""
import hashlib, importlib.metadata as im, inspect, os, sys
import numba, numpy as np
from scipy.sparse import coo_matrix
from umap import umap_

INSTALLED = umap_.smooth_knn_dist


def _build():
    src = inspect.getsource(INSTALLED.py_func)
    src = src[src.index('def smooth_knn_dist'):]
    ns = {k: getattr(umap_, k) for k in ('NPY_INFINITY', 'SMOOTH_K_TOLERANCE', 'MIN_K_DIST_SCALE')}
    ns.update(np=np, numba=numba)
    exec(src, ns)
    fn = numba.njit(locals={'psum': numba.types.float32, 'lo': numba.types.float32, 'mid': numba.types.float32, 'hi': numba.types.float32}, fastmath=False)(ns['smooth_knn_dist'])
    return fn, hashlib.md5(src.encode('utf-8')).hexdigest()


FIXED, SOURCE_MD5 = _build()
umap_.smooth_knn_dist = FIXED

_REVIEW = np.array([[0., 6.9, 7.1, 7.4, 7.9, 8.5, 9.0, 9.6, 10.0, 10.5]], dtype=np.float32)
# (name, distances, k): the first two have a correct sigma above 1, so the fault shows; the third is below 1, where the installed function is right
CASES = (('distances_0_to_14_k15', np.arange(15, dtype=np.float32)[None, :], 15.0),
         ('example_k10', _REVIEW, 10.0),
         ('example_over_20_k10', _REVIEW / 20, 10.0))


def _stop(why):
    sys.exit(f'umap_sigma_fix guard: {why}. This process stops: the integration must not run with a wrong sigma')


def _bindings():
    """every loaded module that holds one of the two names: (those seen, those bound to anything else)"""
    seen, stale = [], []
    for name, m in list(sys.modules.items()):
        d = getattr(m, '__dict__', None)
        if not isinstance(d, dict):
            continue
        s_, f_ = d.get('smooth_knn_dist'), d.get('fuzzy_simplicial_set')
        if s_ is None and f_ is None:
            continue
        seen.append(name)
        if (s_ is not None and s_ is not FIXED) or (f_ is not None and f_ is not umap_.fuzzy_simplicial_set):
            stale.append(name)
    return sorted(seen), sorted(stale)


def guard():
    """Stops the process unless: the override is in place, also in fuzzy_simplicial_set's globals and in every loaded module holding the name; fastmath is off;
    on each known input the compiled value is finite and equals .py_func exactly; and fuzzy_simplicial_set itself returns that value. Returns the record for the run settings."""
    now = umap_.smooth_knn_dist
    if now is not FIXED:
        _stop('umap.umap_.smooth_knn_dist is not the override')
    if umap_.fuzzy_simplicial_set.__globals__.get('smooth_knn_dist') is not FIXED:
        _stop("fuzzy_simplicial_set's globals do not hold the override")
    if now.targetoptions.get('fastmath') is not False:
        _stop(f'the override is compiled with fastmath {now.targetoptions.get("fastmath")!r}')
    seen, stale = _bindings()
    if stale:
        _stop(f'these loaded modules hold another smooth_knn_dist or fuzzy_simplicial_set: {stale}')
    cases = []
    for name, d, k in CASES:
        sig, rho = now(d, k, local_connectivity=1.0)
        psig, prho = now.py_func(d, k, local_connectivity=1.0)
        isig, _ = INSTALLED(d, k, local_connectivity=1.0)
        if not np.isfinite(sig).all():
            _stop(f'{name}: sigma {sig.tolist()} is not finite')
        if not (np.array_equal(sig, psig) and np.array_equal(rho, prho)):
            _stop(f'{name}: compiled sigma {sig.tolist()}, rho {rho.tolist()}; .py_func sigma {psig.tolist()}, rho {prho.tolist()}')
        cases.append(dict(case=name, k=k, sigma_fixed=repr(float(sig[0])), sigma_py_func=repr(float(psig[0])), sigma_installed=repr(float(isig[0]))))
    # the function every caller goes through
    n = 15
    idx = (np.arange(n)[None, :] + np.arange(n)[:, None]) % n
    dst = np.tile(np.arange(n, dtype=np.float32), (n, 1))
    g, s2, _ = umap_.fuzzy_simplicial_set(coo_matrix(([], ([], [])), shape=(n, 1)), n, None, None, knn_indices=idx, knn_dists=dst, set_op_mix_ratio=1.0, local_connectivity=1.0)
    want = now.py_func(dst, float(n), local_connectivity=1.0)[0]
    if not (np.isfinite(s2).all() and np.array_equal(s2, want)):
        _stop(f'fuzzy_simplicial_set returned sigma {s2[:3].tolist()}, .py_func {want[:3].tolist()}')
    me = os.path.abspath(__file__)
    return dict(module=me, module_md5=hashlib.md5(open(me, 'rb').read()).hexdigest(), guard_passed=True, override='umap.umap_.smooth_knn_dist, fastmath=False',
                numba=im.version('numba'), llvmlite=im.version('llvmlite'), umap_learn=im.version('umap-learn'), function_source_md5=SOURCE_MD5,
                installed_fault_present=bool(any(c['sigma_installed'] == 'inf' for c in cases)), cases=cases, route_share_of_weights_exactly_1=round(float((g.data == 1).mean()), 4), modules_holding_the_names=seen)
