#!/usr/bin/env python
# Code release header: converts the two release objects to the CZ CELLxGENE Discover schema 7.1.0 (validator cellxgene-schema 7.1.1).
# Input: all_cells_integrated.h5ad and luminal_cells_integrated.h5ad, opened read-only, gene_id_remapping.csv from a step not included, and the rules of cxg_mappings.py. Outputs are never overwritten.
# It adds the schema fields from per-sample terms and the cell-type map. The var index becomes the Ensembl gene ID.
# A PENDING value of PIPELINE_NEEDED in cxg_mappings.py, or a value of the wrong form, stops the run before an object is read. So does any other failed pre-run check. Methods: Conversion to the CELLxGENE schema.
"""Convert the released atlas AnnData objects (S4: all_cells_integrated.h5ad and
luminal_cells_integrated.h5ad) to the CZ CELLxGENE Discover schema (schema 7.1.0, validator
cellxgene-schema 7.1.1).

Usage (validator environment, which provides anndata 0.11.4 and the pinned ontologies):
  PY=<CXG_ENV>/bin/python
  # small test (two libraries):
  $PY convert_to_cellxgene.py luminal --out-dir <PROJECT_DIR>/final_release_v1/smalltests/S5/<run> --libraries <L1>,<L2>
  # full run (after the small test was checked):
  $PY convert_to_cellxgene.py luminal     # -> cellxgene_schema_v4/luminal_cells_cellxgene_v1.h5ad
  $PY convert_to_cellxgene.py all_cells   # -> cellxgene_schema_v4/all_cells_cellxgene_v1.h5ad

Safety: the S4 objects are opened read-only (h5py mode 'r') and never modified. The new h5ad file is
written element by element (anndata.io.write_elem for obs/var/uns/obsm/varm; chunked CSR streaming
with h5py for X and raw.X) to '<name>.partial' and renamed at the end. Nothing is deleted or
overwritten: an existing output or '.partial' file stops the run; the '.partial' file is
opened in mode 'w-'. Every value that is not decided (cxg_mappings.PENDING, cxg_mappings.PIPELINE_NEEDED)
or has the wrong form (cxg_mappings.check_forms) stops the run before any object is read; the X encoding and
the rules of every S4 column, layer and uns key are checked from the slot keys and column lists before any
value is read. A full run (no --libraries) also stops, before any folder is made, while the author has not confirmed
the 'keep' rules (cxg_mappings.OBS_KEEP_CONFIRMED).

What the conversion does (see cellxgene_field_mapping.csv for details):
  obs  - adds the schema fields (assay, cell type, development stage, disease, donor, ethnicity, sex,
         suspension type, tissue, tissue type, is_primary_data) from the per-library terms and the
         cell-type map (cxg_mappings); the S4 obs columns follow their rules (cell_type and
         cell_type_fine become author_cell_type and author_cell_type_fine; sample-level columns are
         checked equal to sample_metadata.csv); sample-level fields are added from sample_metadata.csv.
  var  - index = Ensembl gene id; gene symbol kept as 'gene_symbol' and the Cell Ranger id as
         'cellranger_gene_id'; 'feature_is_filtered' added. Cell Ranger ids retired in GENCODE v48 (the
         validator reference) are re-labeled with their one-to-one GENCODE v48 successor
         ('gene_id_remapped' = True) or removed from all matrices, as decided in gene_id_remapping.csv.
  raw  - raw.X = the S4 X (raw UMI counts; checked integer and non-negative), float32.
  X    - log1p of counts per 10,000 (natural log) computed from the S4 raw counts, each cell's total
         taken over all Cell Ranger genes before the removal of retired ids, float32.
  layers - not written; each S4 layer needs a rule in cxg_mappings.S4_LAYER_RULES.
  obsm - every S4 embedding of both seeds, unchanged; varm subset to the retained genes.
  obsp - none (the S4 objects have no obsp); an S4 obsp stops the run.
  uns  - title, organism_ontology_term_id, default_embedding, X_approximate_distribution,
         batch_condition, color palettes re-keyed to the kept columns, and the public-text
         'cellxgene_conversion' record.
"""
import argparse
import datetime
import json
import os
import re
import sys
import time

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem, write_elem

sys.dont_write_bytecode = True  # no __pycache__ on a direct run
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cxg_mappings as M  # noqa: E402

# the title is cxg_mappings.TITLE (the article title with the prefix 'All cells - ' or
# 'Luminal epithelial cells - '; check_forms stops unless it is a dict with one str title per dataset)
DEFAULT_EMBEDDING = M.DEFAULT_EMBEDDING
OBSM_KEY_RX = M.OBSM_KEY_RX
CHUNK_ROWS = 20000
MAX_SPAN_ROWS = 60000
SMALLTEST_NEED_BYTES = 2 * 10 ** 9  # [estimate] for two libraries
NEEDED = M.PIPELINE_NEEDED  # every undecided value of the pipeline stops first
# The only text read as a missing value in copied fields; literal 'None', 'NA', 'nan' and 'NaN' are values
# and are kept as they are.
MISSING_TEXT = ("",)
LITERAL_NA_TOKENS = ("None", "NA", "nan", "NaN")  # listed in the log when a written column holds one


def gencode_file():
    import cellxgene_schema
    return os.path.join(os.path.dirname(cellxgene_schema.__file__), "gencode_files", "genes_homo_sapiens.csv.gz")


def log(msg):
    print(f"[{datetime.datetime.now():%Y-%m-%d %H:%M:%S}] {msg}", flush=True)


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def missing_mask(values):
    """True for a real missing value (None, NaN, pd.NA) or '' (MISSING_TEXT); False for every other value,
    including the literal strings 'None', 'NA', 'nan' and 'NaN'."""
    s = pd.Series(values).astype("object")
    return (s.isna() | s.isin(MISSING_TEXT)).values


def as_cat(values, categories=None):
    """Categorical with str categories. Only real missing values and '' become missing; literal strings
    such as 'None', 'NA', 'nan' and 'NaN' are kept as they are."""
    s = pd.Series(values).astype("object")
    miss = missing_mask(s)
    s = s.astype(str).where(~miss, None).astype("object")
    c = pd.Categorical(s, categories=categories) if categories is not None else pd.Categorical(s)
    return c.remove_unused_categories() if categories is None else c


def keep_cat(series):
    """Keep an S4 obs column: categorical columns keep their labels and category order; text columns become
    categorical with the same values; numeric and boolean columns are kept as they are. Only real missing
    values and '' become missing; literal 'None' / 'NA' / 'nan' / 'NaN' stay values."""
    if isinstance(series.dtype, pd.CategoricalDtype):
        miss = missing_mask(series)
        cats = [str(c) for c in series.cat.categories if str(c) not in MISSING_TEXT]
        vals = series.astype(str).where(~miss, None).astype("object")
        return pd.Categorical(vals.values, categories=cats).remove_unused_categories()
    if pd.api.types.is_bool_dtype(series.dtype) or pd.api.types.is_numeric_dtype(series.dtype):
        return series.values
    return as_cat(series.astype("object").values)


def same_column(a, b):
    """Whole-column identity of a source column and its written copy (labels unchanged). Missing means a
    real missing value or ''; literal 'None' / 'NA' / 'nan' / 'NaN' are compared as text."""
    a, b = pd.Series(a), pd.Series(b)
    if pd.api.types.is_numeric_dtype(a.dtype) and pd.api.types.is_numeric_dtype(b.dtype):
        # missing masks first, then the non-missing values (np.array_equal(...,
        # equal_nan=True) raised a TypeError on a nullable boolean column with NA)
        ma, mb = a.isna().to_numpy(), b.isna().to_numpy()
        return len(a) == len(b) and np.array_equal(ma, mb) and bool((a[~ma].to_numpy() == b[~mb].to_numpy()).all())
    ma, mb = missing_mask(a), missing_mask(b)
    if len(a) != len(b) or not np.array_equal(ma, mb):
        return False
    return bool((a.astype(str).values[~ma] == b.astype(str).values[~mb]).all())


def literal_tokens(values):
    """The literal 'None' / 'NA' / 'nan' / 'NaN' strings a written column holds (listed in the log)."""
    s = pd.Series(values).astype("object")
    return sorted(set(s[~missing_mask(s)].astype(str)) & set(LITERAL_NA_TOKENS))


def text_value(v):
    """Exact text of one sample-level value: '' for a real missing value or ''; every other value as str(v),
    not stripped and not compared as a number (for example '01' and '1', or '7.0' and '7', are
    different texts). 'None' / 'NA' / 'nan' / 'NaN' are text."""
    if v is None or (not isinstance(v, str) and pd.isna(v)):
        return ""
    return str(v)


def check_sample_level_column(values, lib, ref, name):
    """STOP unless an S4 obs column named like a sample_metadata.csv column holds one value per library and
    that value is, as exact text, the library's sample_metadata.csv value. values and lib are per cell (same
    index); ref is the sample_metadata.csv column indexed by library."""
    per_lib = pd.Series(values).groupby(lib, observed=True)
    if (per_lib.nunique(dropna=False) > 1).any():
        M.stop(f"sample-level column '{name}' has more than one value within a library")
    first = per_lib.first()
    bad = [x for x in first.index if text_value(first[x]) != text_value(ref.loc[x])]
    if bad:
        M.stop(f"sample-level column '{name}' differs from sample_metadata.csv (exact text) for {bad}")


def read_uns_tolerant(elem, path):
    """read_elem, but entries stored with the 'null' encoding (None, written by a newer anndata
    and unreadable by anndata 0.11) are skipped; returns None for an empty result."""
    if elem.attrs.get("encoding-type") == "null":
        log(f"{path} = None ('null' encoding) not copied")
        return None
    if isinstance(elem, h5py.Group) and elem.attrs.get("encoding-type") == "dict":
        out = {}
        for k in elem.keys():
            v = read_uns_tolerant(elem[k], f"{path}/{k}")
            if v is not None:
                out[k] = v
        if not out:
            log(f"{path} empty after removing 'null' entries; not copied")
            return None
        return out
    return read_elem(elem)


def gencode_ids():
    g = pd.read_csv(gencode_file(), header=None)
    return set(g[0].astype(str).str.split(".").str[0])


def subset_rows(v, rows, n_src):
    if len(rows) == n_src and np.array_equal(rows, np.arange(n_src)):
        return v
    return v.iloc[rows] if isinstance(v, pd.DataFrame) else v[rows]


# ---------------------------------------------------------------------------
# obs
# ---------------------------------------------------------------------------
def build_obs(obs, sm, per):
    for c in (M.S4_LIBRARY_COLUMN, M.S4_CELL_TYPE_COLUMN):
        if c not in obs.columns:
            M.stop(f"source obs has no column '{c}'")
    lib = obs[M.S4_LIBRARY_COLUMN].astype(str)
    missing = set(lib.unique()) - set(sm.index)
    if missing:
        M.stop(f"libraries missing from sample_metadata.csv: {sorted(missing)}")
    smi = sm.loc[sorted(lib.unique())]
    sample_cols = set(sm.columns)

    # every S4 obs column needs a rule (no column is lost or kept silently)
    rules, pending = {}, []
    for c in obs.columns:
        r = M.obs_rule(c, sample_cols)
        if r is None:
            pending.append(c)
        elif r[0] not in M.OBS_ACTIONS:
            M.stop(f"obs rule for '{c}' has unknown action {r[0]}")
        else:
            rules[c] = r
    if pending:
        M.stop(f"source obs columns without a rule (PENDING in cxg_mappings.S4_OBS_RULES): {pending}")

    def by_lib(series_per_lib):
        m = dict(zip(series_per_lib.index, series_per_lib.values))
        return lib.map(m).values

    new = pd.DataFrame(index=obs.index.astype(str))
    new.index.name = None

    # --- schema fields
    new["assay_ontology_term_id"] = as_cat(by_lib(per["assay_ontology_term_id"]))
    ct = obs[M.S4_CELL_TYPE_COLUMN].astype(str)
    unmapped = set(ct.unique()) - set(M.CELL_TYPE_MAP)
    if unmapped:
        M.stop(f"cell type labels without a CL term: {sorted(unmapped)}")
    new["cell_type_ontology_term_id"] = as_cat(ct.map({k: v[0] for k, v in M.CELL_TYPE_MAP.items()}).values)
    new["development_stage_ontology_term_id"] = as_cat(by_lib(per["development_stage_ontology_term_id"]))
    new["disease_ontology_term_id"] = as_cat(by_lib(per["disease_ontology_term_id"]))
    new["donor_id"] = as_cat(by_lib(per["donor_id"]))
    new["is_primary_data"] = np.full(len(new), M.is_primary_data_value(), dtype=bool)
    new["self_reported_ethnicity_ontology_term_id"] = as_cat(by_lib(per["self_reported_ethnicity_ontology_term_id"]))
    new["sex_ontology_term_id"] = as_cat(np.repeat(M.CONSTANT_OBS["sex_ontology_term_id"][0], len(new)))
    new["suspension_type"] = as_cat(np.repeat(M.CONSTANT_OBS["suspension_type"][0], len(new)))
    new["tissue_ontology_term_id"] = as_cat(by_lib(per["tissue_ontology_term_id"]))
    new["tissue_type"] = as_cat(np.repeat(M.CONSTANT_OBS["tissue_type"][0], len(new)))
    assert list(new.columns) == M.SCHEMA_OBS

    # --- S4 obs columns by rule; labels unchanged; sample-level columns equal to sample_metadata.csv
    renamed = {}
    for c in obs.columns:
        action, new_name, _ = rules[c]
        if action == "drop":
            continue
        if new_name != c and new_name in sample_cols:  # never checked against that column
            M.stop(f"source obs column '{c}' would be renamed to '{new_name}', a sample_metadata.csv column name")
        if new_name in new.columns or new_name in M.RESERVED_OBS:
            M.stop(f"source obs column '{c}' would be written as '{new_name}', a schema or reserved name")
        new[new_name] = keep_cat(obs[c])
        renamed[c] = new_name
        if not same_column(obs[c].values, new[new_name].values):
            M.stop(f"values changed in obs column '{c}'")
        # every kept or renamed S4 column named like a sample_metadata.csv column, also one with an
        # S4_OBS_RULES entry, equals sample_metadata.csv per library as exact text
        if c in sample_cols:
            check_sample_level_column(obs[c], lib, smi[c], c)

    # --- sample-level fields added from sample_metadata.csv (authoritative)
    for nc in M.sample_fields_to_add(sm.columns):
        if nc in M.SCHEMA_OBS or nc in M.RESERVED_OBS:  # never skipped silently (build_field_mapping.py, not included, stops too)
            M.stop(f"sample field '{nc}' is a schema or reserved obs name")
        if nc in new.columns:
            continue  # already kept from the S4 obs and checked above (a rename into this name stops above)
        vals = by_lib(smi[nc])
        if nc in M.NUMERIC_SAMPLE_FIELDS:
            try:
                new[nc] = pd.to_numeric(pd.Series(vals).replace("", np.nan)).astype("float64").values
            except (TypeError, ValueError) as e:
                M.stop(f"sample field '{nc}' holds text that is not a number ({e}); it is not turned into a missing value")
        else:
            new[nc] = as_cat(vals)
    for c in new.columns:
        if c not in M.SCHEMA_OBS and literal_tokens(new[c].values):
            log(f"obs column '{c}' keeps the literal values {literal_tokens(new[c].values)} as text (not missing)")
    return new, renamed


def rekey_colors(uns_src, obs, new_obs, renamed):
    """'<column>_colors' palettes follow their column (kept, renamed or dropped), same category order."""
    out = {}
    for k, v in uns_src.items():
        if not k.endswith("_colors"):
            continue
        old_col = k[: -len("_colors")]
        if old_col not in renamed:
            log(f"color key {k} dropped (column not kept)")
            continue
        new_col = renamed[old_col]
        new_key = f"{new_col}_colors"
        if not isinstance(obs[old_col].dtype, pd.CategoricalDtype) or not isinstance(new_obs[new_col].dtype, pd.CategoricalDtype):
            log(f"color key {k} dropped (column not categorical)")
            continue
        cmap = dict(zip([str(c) for c in obs[old_col].cat.categories], list(v)))
        new_cats = [str(c) for c in new_obs[new_col].cat.categories]
        if all(c in cmap for c in new_cats):
            out[new_key] = np.array([cmap[c] for c in new_cats], dtype=object)
            log(f"colors {k} -> {new_key} ({len(new_cats)} categories)")
        else:
            log(f"color key {k} dropped (categories changed)")
    return out


# ---------------------------------------------------------------------------
# matrices
# ---------------------------------------------------------------------------
def stream_csr(src_grp, dst_file, dst_path, rows, col_remap, n_cols_out, dtype, mode):
    """Copy selected rows (sorted int array of source row ids) of a source CSR group into a new
    CSR group, keeping only columns with col_remap >= 0 (remap preserves column order).
    mode 'raw': values copied, and counted when not integer or negative;
    mode 'lognorm': log1p(value / row total * 10,000), the row total over all source columns."""
    if src_grp.attrs.get("encoding-type") != "csr_matrix":
        M.stop(f"source matrix is {src_grp.attrs.get('encoding-type')}, not csr_matrix")
    if mode not in ("raw", "lognorm"):
        raise ValueError(mode)
    ip = src_grp["indptr"][...].astype(np.int64)
    sdata, sind = src_grp["data"], src_grp["indices"]
    g = dst_file.create_group(dst_path)
    d_data = g.create_dataset("data", shape=(0,), maxshape=(None,), dtype=dtype, chunks=(1 << 20,))
    d_ind = g.create_dataset("indices", shape=(0,), maxshape=(None,), dtype=np.int32, chunks=(1 << 20,))
    indptr = np.zeros(len(rows) + 1, dtype=np.int64)
    gene_sum = np.zeros(n_cols_out, dtype=np.float64)
    gene_nnz = np.zeros(n_cols_out, dtype=np.int64)
    stats = dict(rows_all_zero=0, non_integer=0, negative=0, explicit_zero=0, nnz=0)
    pos = 0
    t0 = time.time()
    c0 = 0
    n_chunks = 0
    while c0 < len(rows):
        c1 = min(c0 + CHUNK_ROWS, len(rows))
        while c1 - c0 > 1 and rows[c1 - 1] - rows[c0] >= MAX_SPAN_ROWS:  # bound the source span read at once
            c1 = c0 + (c1 - c0) // 2
        sel = rows[c0:c1]
        r0, r1 = int(sel[0]), int(sel[-1]) + 1
        a, b = ip[r0], ip[r1]
        data = sdata[a:b]
        ind = sind[a:b]
        lens = ip[sel + 1] - ip[sel]
        starts = ip[sel] - a
        if len(sel) == r1 - r0:  # contiguous block
            take = None
        else:
            offs = np.repeat(starts - np.concatenate([[0], np.cumsum(lens)[:-1]]), lens)
            take = offs + np.arange(lens.sum())
        if take is not None:
            data = data[take]
            ind = ind[take]
        data = data.astype(np.float64)
        if mode == "raw":
            stats["non_integer"] += int((data != np.round(data)).sum())
            stats["negative"] += int((data < 0).sum())
            stats["explicit_zero"] += int((data == 0).sum())
        else:
            row_of = np.repeat(np.arange(len(sel)), lens)
            tot = np.bincount(row_of, weights=data, minlength=len(sel))[row_of]
            data = np.log1p(np.divide(data, tot, out=np.zeros_like(data), where=tot > 0) * 1e4)
        newc = col_remap[ind]
        keep = newc >= 0
        data = data[keep].astype(dtype)
        newc = newc[keep].astype(np.int32)
        cm = np.concatenate([[0], np.cumsum(keep, dtype=np.int64)])
        bounds = np.concatenate([[0], np.cumsum(lens)])
        counts = cm[bounds[1:]] - cm[bounds[:-1]]
        indptr[c0 + 1:c0 + 1 + len(sel)] = pos + np.cumsum(counts)
        n = len(data)
        d_data.resize((pos + n,))
        d_ind.resize((pos + n,))
        d_data[pos:pos + n] = data
        d_ind[pos:pos + n] = newc
        pos += n
        stats["rows_all_zero"] += int((counts == 0).sum())
        if mode == "raw":
            gene_sum += np.bincount(newc, weights=data, minlength=n_cols_out)
            gene_nnz += np.bincount(newc, minlength=n_cols_out)
        n_chunks += 1
        if n_chunks % 10 == 1:
            log(f"  {dst_path}: {c0 + len(sel):,}/{len(rows):,} rows ({time.time() - t0:.0f}s)")
        c0 = c1
    g.create_dataset("indptr", data=indptr)
    g.attrs["encoding-type"] = "csr_matrix"
    g.attrs["encoding-version"] = "0.1.0"
    g.attrs["shape"] = np.array([len(rows), n_cols_out], dtype=np.int64)
    stats["nnz"] = int(pos)
    return stats, gene_sum, gene_nnz


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main(dataset, out_dir, libraries_arg):
    if dataset not in M.OUT_NAME:
        M.stop(f"dataset must be one of {list(M.OUT_NAME)}")
    M.require_decided(NEEDED)
    M.check_forms(NEEDED)  # the forms of the needed values
    M.check_uns_rule_names()  # kept S4 uns keys against the schema's uns names
    versions = M.check_versions()
    smalltest = bool(libraries_arg)
    if not smalltest:  # a full run waits for the author's confirmation of the 'keep' rules
        M.check_full_run_confirmed()
    M.check_out_dir(out_dir, smalltest)
    src_path = M.PATHS[dataset]
    out_path = os.path.join(out_dir, M.OUT_NAME[dataset])
    tmp_path = out_path + ".partial"
    for p in (out_path, tmp_path):
        if os.path.exists(p):
            M.stop(f"{p} exists; nothing is overwritten or removed by this script")
    M.check_free_space(out_dir, SMALLTEST_NEED_BYTES if smalltest else M.EXPECTED_OUTPUT_BYTES[dataset])
    st0 = os.stat(src_path)
    log(f"dataset={dataset} source={src_path} (read-only, size {st0.st_size:,}) -> {out_path}; versions {versions}")

    from cellxgene_ontology_guide.ontology_parser import OntologyParser
    P = OntologyParser(schema_version=f"v{M.SCHEMA_VERSION}")

    sm = M.read_sample_metadata()
    libs = M.parse_libraries(libraries_arg, set(sm.index))
    per = M.per_sample_schema_fields(sm, P.get_term_label_to_id_map("HsapDv"))

    # ---------------- slot keys, column lists and the X encoding first (attributes only) ----------------
    with h5py.File(src_path, "r") as f:
        if "raw" in f:
            M.stop("the source object has a raw slot; raw counts are expected in X")
        if "obsp" in f and len(f["obsp"].keys()):
            M.stop(f"the source object has obsp {list(f['obsp'].keys())}; none is expected")
        # a CSC or dense X stops here, before the '.partial' file exists
        x_enc = f["X"].attrs.get("encoding-type") if "X" in f else "missing"
        if x_enc != "csr_matrix":
            M.stop(f"source X is {x_enc}, not csr_matrix")
        layer_keys = list(f["layers"].keys()) if "layers" in f else []
        uns_keys = list(f["uns"].keys()) if "uns" in f else []
        obs_cols = [str(c) for c in f["obs"].attrs["column-order"]]
        var_cols = [str(c) for c in f["var"].attrs["column-order"]]

    # every S4 obs / var column, layer and uns key needs a rule with a known action; checked before any value
    # is read, so an unreadable entry with a 'drop' rule is never read
    sample_cols = set(sm.columns)
    pend, bad_action = [], []
    for c in obs_cols:
        r = M.obs_rule(c, sample_cols)
        if r is None:
            pend.append(f"obs/{c}")
        elif r[0] not in M.OBS_ACTIONS:
            bad_action.append(f"obs/{c}: {r[0]}")
    for c in var_cols:
        if c in (M.S4_GENE_ID_SOURCE, M.S4_GENE_SYMBOL_SOURCE):
            continue
        r = M.S4_VAR_RULES.get(c)
        if r is None:
            pend.append(f"var/{c}")
        elif r[0] not in M.VAR_ACTIONS:
            bad_action.append(f"var/{c}: {r[0]}")
    for k in layer_keys:
        r = M.S4_LAYER_RULES.get(k)
        if r is None:
            pend.append(f"layers/{k}")
        elif r[0] not in M.LAYER_ACTIONS:
            bad_action.append(f"layers/{k}: {r[0]} (layers are not written; only 'drop')")
    for k in uns_keys:
        if k.endswith("_colors"):
            continue
        r = M.S4_UNS_RULES.get(k)
        if r is None:
            pend.append(f"uns/{k}")
        elif r[0] not in M.UNS_ACTIONS:
            bad_action.append(f"uns/{k}: {r[0]}")
    if pend:
        M.stop(f"source columns / layers / uns keys without a rule (PENDING in cxg_mappings): {pend}")
    if bad_action:
        M.stop(f"obs, var, layer or uns rules with an unknown action: {bad_action}")
    uns_read = [k for k in uns_keys if k.endswith("_colors") or M.S4_UNS_RULES[k][0] == "keep"]

    # ---------------- read source metadata (read-only) ----------------
    with h5py.File(src_path, "r") as f:
        obs = read_elem(f["obs"])
        var = read_elem(f["var"])
        obsm_src = {k: read_elem(f["obsm"][k]) for k in f["obsm"].keys()} if "obsm" in f else {}
        varm_src = {k: read_elem(f["varm"][k]) for k in f["varm"].keys()} if "varm" in f else {}
        uns_src = {k: read_uns_tolerant(f["uns"][k], f"uns/{k}") for k in uns_read}
        uns_src = {k: v for k, v in uns_src.items() if v is not None}
        src_shape = tuple(int(x) for x in f["X"].attrs["shape"])
        x_dtype = f["X"]["data"].dtype
    if [str(c) for c in obs.columns] != obs_cols or [str(c) for c in var.columns] != var_cols:
        M.stop("obs / var columns read differ from the column lists checked above")
    log(f"source: {src_shape[0]:,} cells x {src_shape[1]:,} genes; X dtype {x_dtype}; obs cols {obs.shape[1]}; "
        f"var cols {var.shape[1]}; obsm {list(obsm_src)}; layers {layer_keys} (dropped); uns {uns_keys} (read: {uns_read})")
    if not smalltest and M.N_CELLS[dataset] is not None and src_shape[0] != M.N_CELLS[dataset]:
        M.stop(f"{src_shape[0]:,} cells (expected {M.N_CELLS[dataset]:,})")

    # ---------------- var ----------------
    ids, sym = M.s4_gene_ids_symbols(var)
    if src_shape[1] != len(ids):
        M.stop(f"X has {src_shape[1]} columns but var has {len(ids)} genes")
    if dataset == "luminal":
        with h5py.File(M.PATHS["all_cells"], "r") as fa:
            a_ids, a_sym = M.s4_gene_ids_symbols(read_elem(fa["var"]))
        if not (np.array_equal(ids, a_ids) and np.array_equal(sym, a_sym)):
            M.stop("gene order or names differ between the Luminal and all-cells objects")
    v48 = gencode_ids()
    remap = pd.read_csv(os.path.join(out_dir, "gene_id_remapping.csv"), dtype=str, keep_default_na=False).set_index("cellranger_gene_id")
    retired = [i for i in ids if i not in v48]
    if set(retired) != set(remap.index):
        M.stop("gene_id_remapping.csv does not match the retired ids of this object; run build_gene_id_remap.py")
    to_new = {i: i for i in ids if i in v48}
    to_new.update({i: remap.loc[i, "candidate_id"] for i in retired if remap.loc[i, "decision"] == "remap"})
    ok = np.array([i in to_new for i in ids])
    new_ids = np.array([to_new[i] for i in ids[ok]])
    if not (pd.Index(new_ids).is_unique and set(new_ids) <= v48):
        M.stop("remapped gene ids are not unique or not in GENCODE v48")
    col_remap = np.full(len(ids), -1, dtype=np.int64)
    col_remap[ok] = np.arange(ok.sum())
    n_remapped = int(sum(1 for i in ids[ok] if i not in v48))
    log(f"genes kept: {ok.sum():,} (of which {n_remapped} retired identifiers relabeled with their GENCODE v48 successor); removed: {(~ok).sum():,}")
    new_var = pd.DataFrame(index=pd.Index(new_ids, name=None))
    new_var["gene_symbol"] = np.asarray(sym[ok], dtype=str).astype(object)
    new_var["cellranger_gene_id"] = ids[ok].astype(object)
    new_var["gene_id_remapped"] = np.array([i not in v48 for i in ids[ok]], dtype=bool)
    new_var["feature_is_filtered"] = np.zeros(ok.sum(), dtype=bool)
    pend = []
    for c in var.columns:
        if c in (M.S4_GENE_ID_SOURCE, M.S4_GENE_SYMBOL_SOURCE):
            continue
        rule = M.S4_VAR_RULES.get(c)
        if rule is None:
            pend.append(c)
            continue
        action, new_name, _ = rule
        if action == "drop":
            if c in M.S4_VAR_DROP_IF_CONSTANT and set(var[c].astype(str)) != {M.S4_VAR_DROP_IF_CONSTANT[c]}:
                M.stop(f"var column '{c}' is not the constant {M.S4_VAR_DROP_IF_CONSTANT[c]!r}")
            continue
        if new_name in new_var.columns or new_name in M.RESERVED_VAR:
            M.stop(f"var column '{c}' would be written as '{new_name}', a used or reserved name")
        new_var[new_name] = var[c].values[ok]
    if pend:
        M.stop(f"source var columns without a rule (PENDING in cxg_mappings.S4_VAR_RULES): {pend}")
    raw_var = pd.DataFrame(index=new_var.index.copy())

    # ---------------- rows (small test: the named libraries only) ----------------
    lib_all = obs[M.S4_LIBRARY_COLUMN].astype(str) if M.S4_LIBRARY_COLUMN in obs.columns else None
    if lib_all is None:
        M.stop(f"source obs has no column '{M.S4_LIBRARY_COLUMN}'")
    rows = np.arange(src_shape[0], dtype=np.int64) if libs is None else np.flatnonzero(lib_all.isin(libs).values).astype(np.int64)
    if len(rows) == 0:
        M.stop(f"no cells of {libs} in this object")
    obs_sel = obs.iloc[rows]

    # ---------------- obs ----------------
    new_obs, renamed = build_obs(obs_sel, sm, per)
    log(f"obs: {len(new_obs):,} cells, {new_obs.shape[1]} columns")

    # ---------------- obsm / varm ----------------
    bad = [k for k in obsm_src if not re.match(OBSM_KEY_RX, k)]
    if bad:
        M.stop(f"obsm keys not in the schema form X_<suffix>: {bad}")
    if DEFAULT_EMBEDDING[dataset] not in obsm_src:
        M.stop(f"default embedding {DEFAULT_EMBEDDING[dataset]} not in obsm {list(obsm_src)}")
    obsm = {k: subset_rows(v, rows, src_shape[0]) for k, v in obsm_src.items()}
    varm = {k: (v.iloc[ok] if isinstance(v, pd.DataFrame) else v[ok]) for k, v in varm_src.items()}

    # ---------------- uns ----------------
    uns = {
        "title": M.TITLE[dataset],
        "organism_ontology_term_id": M.ORGANISM_ONTOLOGY_TERM_ID,
        "default_embedding": DEFAULT_EMBEDDING[dataset],
        "X_approximate_distribution": "normal",
        "batch_condition": np.array(["study_id"], dtype=object),
    }
    if "study_id" not in new_obs.columns:
        M.stop("batch_condition names study_id, which is not in obs")
    for k in uns_keys:
        if k in M.S4_UNS_RULES and M.S4_UNS_RULES[k][0] == "keep":
            if k in uns or k not in uns_src:
                M.stop(f"source uns key '{k}' collides with a schema uns key or could not be read")
            uns[k] = uns_src[k]
            if k == "embeddings":  # '(cell_type)' -> '(author_cell_type)' in the all-cell training_labels
                uns[k] = M.rename_training_label_tokens(uns[k], dataset)
    uns.update(rekey_colors(uns_src, obs_sel, new_obs, renamed))
    uns["cellxgene_conversion"] = {
        "description": (f"Converted to the CZ CELLxGENE Discover schema {M.SCHEMA_VERSION} (validator cellxgene-schema "
                        f"v{M.PINNED_PACKAGES['cellxgene-schema']}). X: log(1 + counts per 10,000), natural log. Each cell's "
                        f"total is over all {M.N_RELEASE_GENES:,} Cell Ranger genes, before identifiers retired in GENCODE v48 "
                        "are removed. raw.X: raw UMI counts (Cell Ranger v7.1.0, GRCh38-2020-A). In four studies, "
                        "raw.X excludes barcode-swapped molecules (Methods)."),
        "source_file": M.PUBLIC_SOURCE_NAME[dataset],
        "date": datetime.date.today().isoformat(),
        "sample_metadata": "sample_metadata.csv of the same deposit",
    }
    hits = {k: M.public_text_problems(v) for k, v in list(uns["cellxgene_conversion"].items()) + [("title", uns["title"])]}
    hits = {k: v for k, v in hits.items() if v}
    if hits:
        M.stop(f"public-text scan of uns failed: {hits}")

    # ---------------- write ----------------
    n_out = int(ok.sum())
    with h5py.File(tmp_path, "w-") as fo, h5py.File(src_path, "r") as fs:  # 'w-': never an existing file
        fo.attrs["encoding-type"] = "anndata"
        fo.attrs["encoding-version"] = "0.1.0"
        write_elem(fo, "obs", new_obs)
        write_elem(fo, "var", new_var)
        log("writing raw.X (raw UMI counts = source X)")
        rstats, gsum, gnnz = stream_csr(fs["X"], fo, "raw/X", rows, col_remap, n_out, np.float32, mode="raw")
        log(f"raw.X done: {rstats}")
        if rstats["non_integer"] or rstats["negative"]:
            M.stop(f"source X is not raw counts ({rstats}); {tmp_path} is left in place and not removed by this script")
        write_elem(fo["raw"], "var", raw_var)
        write_elem(fo["raw"], "varm", {})
        fo["raw"].attrs["encoding-type"] = "raw"
        fo["raw"].attrs["encoding-version"] = "0.1.0"
        log("writing X (log1p of counts per 10,000, float32)")
        xstats, _, _ = stream_csr(fs["X"], fo, "X", rows, col_remap, n_out, np.float32, mode="lognorm")
        log(f"X done: {xstats}")
        write_elem(fo, "obsm", obsm)
        write_elem(fo, "varm", varm)
        write_elem(fo, "obsp", {})
        write_elem(fo, "varp", {})
        write_elem(fo, "layers", {})
        write_elem(fo, "uns", uns)

    st1 = os.stat(src_path)
    if (st1.st_size, st1.st_mtime) != (st0.st_size, st0.st_mtime):
        M.stop(f"the source file changed during the run ({src_path})")
    summary = dict(dataset=dataset, smalltest_libraries=libs, n_obs=int(len(new_obs)), n_vars=n_out,
                   genes_removed=int((~ok).sum()), genes_remapped=n_remapped, raw=rstats, X=xstats,
                   obs_columns=list(new_obs.columns), var_columns=list(new_var.columns), uns_keys=list(uns.keys()),
                   obsm_keys=list(obsm.keys()), varm_keys=list(varm.keys()), dropped_layers=layer_keys,
                   versions=versions, source=src_path, source_size=st0.st_size)
    summ_path = os.path.join(out_dir, "logs", f"convert_{dataset}_summary.json")
    if os.path.exists(summ_path):
        M.stop(f"{summ_path} exists")
    with open(summ_path, "x") as fh:  # 'x': never an existing file
        json.dump(summary, fh, indent=1, default=str)
    if os.path.exists(out_path):
        M.stop(f"{out_path} appeared during the run; {tmp_path} kept")
    os.rename(tmp_path, out_path)
    log(f"finished: {out_path} ({os.path.getsize(out_path) / 1e9:.1f} GB)")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("dataset", choices=sorted(M.OUT_NAME))
    ap.add_argument("--out-dir", default=M.PATHS["out_dir"])
    ap.add_argument("--libraries", default="", help="small test: comma-separated library ids (no default)")
    a = ap.parse_args()
    main(a.dataset, a.out_dir, a.libraries)
