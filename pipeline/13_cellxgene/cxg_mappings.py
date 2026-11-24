# Code release header: mapping rules and checks for the CELLxGENE conversion, imported by convert_to_cellxgene.py in this folder and by scripts not included.
# Inputs: sample_metadata.csv (md5 pinned) and, through the converter, the two release objects. Methods: Conversion to the CELLxGENE schema.
# Key settings: CELLxGENE schema 7.1.0 with cellxgene-schema 7.1.1, and the ontology term maps of cell type, disease, tissue, assay, age and ethnicity.
# A constant set to PENDING would stop the run, and none is PENDING. Two values gate only a final check of the output tables, and one gates a full run.
"""Mapping rules used to convert the released atlas AnnData objects (the S4 release objects
all_cells_integrated.h5ad and luminal_cells_integrated.h5ad) to the CZ CELLxGENE Discover schema
(schema 7.1.0, as in the article, Table 2; validator cellxgene-schema 7.1.1).

This module is the single source of truth for
  * the pinned versions, input paths and output folders,
  * ontology mappings (cell type, disease, tissue, assay, development stage, ethnicity),
  * the obs / var / layers / uns plan for the S4 objects (keep / rename / drop / derived),
  * the decided per-library terms and the checks against them,
  * the public-text scan (American spelling, no internal names) of every text written into the
    mapping tables and into uns['cellxgene_conversion'].
`build_mapping_tables.py` (not included) writes the human-readable CSV tables from it, and
`convert_to_cellxgene.py` applies it to the h5ad files.

It holds the decided conversion values, with the S4 layout values written from the S4 objects. The 'keep' rules of
7 obs columns and of the uns key 'embeddings' block a full run until the author confirms them (OBS_KEEP_CONFIRMED).
A value that is still not decided is the constant PENDING. No value of the pipeline list is PENDING. Step 0 of the
run script (not included) (preflight()) and every pipeline script check the whole list PIPELINE_NEEDED first and stop, before any
folder is made, any object is read or any X pass starts, while one of them is PENDING (require_decided)
or a value has the wrong form (check_forms). Two values wait for later steps and block only the final
hand-over of the outputs (check_handover(), HANDOVER_NEEDED): the author's confirmation of the flagged terms
and the reviewed S4 Monocyte count.

Every ontology term id below was checked against the ontology releases pinned by
cellxgene-schema 7.1.1 (cellxgene-ontology-guide 1.11.1): CL 2026-03-26, EFO 3.90.0,
MONDO 2026-05-05, UBERON 2026-04-01, HsapDv 2025-01-23, HANCESTRO 2025-10-14, PATO 2025-05-14.
`build_mapping_tables.py` (not included) re-checks them (label lookup, not deprecated) every time it runs.
"""
import hashlib
import os
import re

import numpy as np
import pandas as pd

PENDING = "PENDING"  # a value that is not decided yet; require_decided() stops on it


def stop(msg, code=2):
    """Stop the run with a message (no fallback, no automatic switch)."""
    print(f"STOP: {msg}", flush=True)
    raise SystemExit(code)


# ----------------------------------------------------------------------------
# Paths (code constants). All outputs under <PROJECT_DIR>.
# ----------------------------------------------------------------------------
REV = "<PROJECT_DIR>"
PATHS = dict(
    # S4 release objects: one object per set
    all_cells=f"{REV}/final_release_objects_v1e/all_cells_integrated.h5ad",
    luminal=f"{REV}/final_release_objects_v1e/luminal_cells_integrated.h5ad",
    # metadata v35d, the version pinned here: 163 rows, no included / exclusion_reason columns,
    # 60 columns (N_SAMPLE_COLUMNS) and the library key library_id
    sample_metadata=f"{REV}/repository_files_v35d_2026-11-22/sample_metadata.csv",
    metadata_dictionary=f"{REV}/repository_files_v35d_2026-11-22/metadata_dictionary.csv",
    # full-run outputs and the S5 small test
    out_dir=f"{REV}/cellxgene_schema_v4",
    smalltest_dir=f"{REV}/final_release_v1/smalltests/S5",
    # first-version inputs, read only: cached Ensembl archive result and GENCODE v32 GTF (no
    # Ensembl query; S5 recomputes the UMI shares only), and the earlier conversion's tables that these results are compared with
    reference_dir=f"{REV}/cellxgene_schema_v2/reference",
    v2_gene_id_remapping=f"{REV}/cellxgene_schema_v2/gene_id_remapping.csv",
    v2_sample_level_fields=f"{REV}/cellxgene_schema_v2/sample_level_schema_fields.csv",
)
# output names as the first-version files (same names, other folder)
OUT_NAME = {"all_cells": "all_cells_cellxgene_v1.h5ad", "luminal": "luminal_cells_cellxgene_v1.h5ad"}
# the deposit names of the S4 objects; written to uns['cellxgene_conversion']['source_file']
PUBLIC_SOURCE_NAME = {"all_cells": "all_cells_integrated.h5ad", "luminal": "luminal_cells_integrated.h5ad"}
# uns['default_embedding']: all cells X_umap_scanvi_seed0 (primary);
# Luminal X_umap_scanorama_seed0; obsm keys as in the S4 objects
DEFAULT_EMBEDDING = {"all_cells": "X_umap_scanvi_seed0", "luminal": "X_umap_scanorama_seed0"}
OBSM_KEY_RX = r"^X_[a-zA-Z][a-zA-Z0-9_.-]*$"  # as cellxgene-schema 7.1.1 validate.py l.1528-1547

# md5 of the sample_metadata.csv the run reads (v35d), computed in the release preparation
# and equal to the md5 recorded for v35d. The three term maps were
# re-checked against this file (ASSAY_MAP, TISSUE_SITE_MAP, ETHNICITY_MAP below). A changed file stops the run.
SAMPLE_METADATA_MD5 = "b8b98492e2c2dd2b03c734c10c1da9e9"

# uns['title'] of the two files (the author's decision: "the submitted article title, prefixed with
# 'All cells - ' or 'Luminal epithelial cells - '"). The article title was taken from the current manuscript and
# equals the submitted title (the dataset title follows the article title).
ARTICLE_TITLE = ("A harmonized single-cell RNA-seq atlas of human localized and metastatic prostate cancers and "
                 "benign tissues")
TITLE_PREFIX = {"all_cells": "All cells - ", "luminal": "Luminal epithelial cells - "}
TITLE = {ds: prefix + ARTICLE_TITLE for ds, prefix in TITLE_PREFIX.items()}

# ----------------------------------------------------------------------------
# Pinned versions (schema 7.1.0 with cellxgene-schema 7.1.1)
# ----------------------------------------------------------------------------
SCHEMA_VERSION = "7.1.0"
PINNED_PACKAGES = {"cellxgene-schema": "7.1.1", "cellxgene-ontology-guide": "1.11.1", "anndata": "0.11.4"}


def check_versions():
    """Stop unless the validator environment has exactly the pinned packages."""
    from importlib.metadata import PackageNotFoundError, version
    bad = {}
    for pkg, want in PINNED_PACKAGES.items():
        try:
            got = version(pkg)
        except PackageNotFoundError:
            got = "not installed"
        if got != want:
            bad[pkg] = (got, want)
    if bad:
        stop(f"package versions differ from the pins (found, pinned): {bad}")
    # pandas is recorded, not pinned (the public-text scan depends on its dtypes)
    return dict(PINNED_PACKAGES, schema=SCHEMA_VERSION, pandas_recorded=pd.__version__)


# ----------------------------------------------------------------------------
# Run guards
# ----------------------------------------------------------------------------
MIN_FREE_BYTES = 200 * 2 ** 30  # free space kept on the disk of <PROJECT_DIR> (as the run wrapper's start guard)
# expected output sizes: first-version files plus both seeds' obsm (about 1.2 / 0.4 GB) [estimate]
EXPECTED_OUTPUT_BYTES = {"all_cells": 24_467_284_169 + 1_200_000_000, "luminal": 8_036_051_793 + 400_000_000}
N_LIBRARIES = 163          # sample_metadata.csv rows (v35d)
N_RELEASE_GENES = 36601    # the 36,601 release gene ids
N_CELLS = {"all_cells": 740975, "luminal": 147667}  # the cell counts of the two S4 objects
N_RETIRED_IDS_V2 = 1142    # from the earlier gene_id_remapping.csv
N_MT_GENES = 13            # the 13 MT- genes


def _is_pending(v):
    if isinstance(v, str):
        return v == PENDING
    if isinstance(v, dict):
        return any(_is_pending(x) for x in v.values())
    if isinstance(v, (list, tuple)):
        return any(_is_pending(x) for x in v)
    return False


def require_decided(names):
    """Stop before any object is read while one of the named module constants is PENDING."""
    pend = []
    for n in names:
        v = eval(n, globals())
        if _is_pending(v):  # for a map, name the PENDING entries
            pend.append(f"{n}{[k for k, x in v.items() if _is_pending(x)]}" if isinstance(v, dict) else n)
    if pend:
        stop("undecided constants (PENDING): " + ", ".join(pend) + " - see the comment at each constant in cxg_mappings.py")


# Every value of the pipeline that was undecided at the start of the release preparation. Each pipeline script calls
# require_decided(PIPELINE_NEEDED) as the first statement of main() and check_forms(PIPELINE_NEEDED) next, and
# step 0 of the run script (not included) calls preflight() before it makes any folder, so no object is read and no X pass starts
# while any of these is PENDING (a map counts as PENDING while one of its entries is) or of the wrong form.
# SAMPLE_METADATA_MD5, CONSTANT_OBS, SAMPLE_FIELDS_ADD, TITLE and the three term maps are decided
# and stay on the list for their form check. The S4 layout values (gene sources, the four *_RULES_FROM
# sentinels, S4_ALL_CELLS_HVG_FLAGS) come from the S4 objects and stay listed for their form check.
S4_LAYOUT_NEEDED = ["S4_OBS_RULES_FROM", "S4_VAR_RULES_FROM", "S4_LAYER_RULES_FROM", "S4_UNS_RULES_FROM"]
TERM_MAPS_NEEDED = ["ASSAY_MAP", "TISSUE_SITE_MAP", "ETHNICITY_MAP"]
PIPELINE_NEEDED = ["SAMPLE_METADATA_MD5", "CONSTANT_OBS", "SAMPLE_FIELDS_ADD", "S4_GENE_ID_SOURCE",
                   "S4_GENE_SYMBOL_SOURCE", "S4_ALL_CELLS_HVG_FLAGS", "TITLE"] + S4_LAYOUT_NEEDED + TERM_MAPS_NEEDED
# Values that do not stop the pipeline but block the final hand-over of the outputs (check_handover()):
# the author's confirmation of the flagged terms and the reviewed S4 Monocyte count.
HANDOVER_NEEDED = ["AUTHOR_REVIEW_CONFIRMED", "MONOCYTE_S4_COUNT_FROM"]
# The value that blocks only a full run: the author's confirmation of the 'keep' rules. check_full_run_confirmed()
# checks it in step 0 of a full run of the run script (not included) and in convert_to_cellxgene.py without
# --libraries. A small test does not check it.
FULL_RUN_NEEDED = ["OBS_KEEP_CONFIRMED"]


def _str_tuples(m, n):
    """Every value of the map m is a tuple of n str values."""
    return isinstance(m, dict) and all(isinstance(x, tuple) and len(x) == n and all(isinstance(y, str) for y in x)
                                       for x in m.values())


# The form of every decided value of PIPELINE_NEEDED. A value of the wrong form stops
# here with its name, before any folder is made or any object is read, instead of failing with a traceback
# later (a TITLE written as a str failed only after both objects were read). check_forms() runs after
# require_decided(), so no value checked here is PENDING.
DECIDED_FORMS = {
    "SAMPLE_METADATA_MD5": (lambda v: isinstance(v, str) and re.fullmatch(r"[0-9a-f]{32}", v) is not None,
                            "a lower-case md5 hex digest (str)"),
    "CONSTANT_OBS": (lambda v: isinstance(v["is_primary_data"], tuple) and len(v["is_primary_data"]) == 2
                     and isinstance(v["is_primary_data"][0], bool) and isinstance(v["is_primary_data"][1], str),
                     "a dict whose 'is_primary_data' is (True or False, note) with a Python bool"),
    "SAMPLE_FIELDS_ADD": (lambda v: (isinstance(v, str) and v == "all")
                          or (isinstance(v, list) and all(isinstance(x, str) for x in v)),
                          "'all' or a list of sample_metadata.csv column names (str)"),
    "S4_GENE_ID_SOURCE": (lambda v: isinstance(v, str) and v != "", "'index' or a var column name (str)"),
    "S4_GENE_SYMBOL_SOURCE": (lambda v: isinstance(v, str) and v != "", "'index' or a var column name (str)"),
    "S4_ALL_CELLS_HVG_FLAGS": (lambda v: isinstance(v, list) and len(set(v)) == len(v)
                               and all(isinstance(x, str) and x != "" for x in v),
                               "a list (possibly empty) of distinct S4 all-cells var column names (str)"),
    "TITLE": (lambda v: isinstance(v, dict) and set(v) == set(OUT_NAME)
              and all(isinstance(x, str) and x.strip() != "" for x in v.values()),
              f"a dict with the keys {sorted(OUT_NAME)} and a non-empty str title for each"),
    **{n: (lambda v: isinstance(v, str) and v.strip() != "", "a non-empty str ('<record>; <md5>')")
       for n in S4_LAYOUT_NEEDED + HANDOVER_NEEDED + FULL_RUN_NEEDED},
    "ASSAY_MAP": (lambda v: _str_tuples(v, 2), "a dict of (term id, note) str pairs"),
    "TISSUE_SITE_MAP": (lambda v: _str_tuples(v, 3), "a dict of (term id, note, review flag) str triples"),
    "ETHNICITY_MAP": (lambda v: _str_tuples(v, 2), "a dict of (term id, note) str pairs"),
}


def check_forms(names):
    """Stop while one of the named decided constants has the wrong form (run after require_decided)."""
    bad = []
    for n in names:
        if n not in DECIDED_FORMS:
            bad.append(f"{n} has no form rule")
            continue
        ok, want = DECIDED_FORMS[n]
        try:
            good = bool(ok(eval(n, globals())))
        except (AttributeError, KeyError, TypeError, ValueError):
            good = False
        if not good:
            bad.append(f"{n} must be {want}")
    if bad:
        stop("decided constants of the wrong form: " + "; ".join(bad) + " - see the comment at each constant in cxg_mappings.py")


def preflight():
    """Step 0 of the run script (not included): checks the constants and package versions only (no file, folder or object)."""
    require_decided(PIPELINE_NEEDED)
    check_forms(PIPELINE_NEEDED)
    check_uns_rule_names()  # kept S4 uns keys against the schema's uns names
    check_review_flags()  # review-flag consistency; PENDING hand-over values do not stop here
    v = check_versions()
    print(f"preflight passed: no PENDING value, every decided value of the expected form; versions {v}", flush=True)


def check_out_dir(out_dir, smalltest):
    """Outputs only under the REV folder; a small test never writes into the full-run folder."""
    real = os.path.realpath(out_dir)
    rev = os.path.realpath(REV)
    if not real.startswith(rev + os.sep):
        stop(f"output folder {out_dir} is outside the allowed project tree (REV in cxg_mappings.py)")
    # a small test writes into a run folder below smalltests/S5/, never into smalltests/S5/ itself (the run
    # name '.' once wrote straight into it)
    if smalltest and not real.startswith(os.path.realpath(PATHS["smalltest_dir"]) + os.sep):
        stop(f"a small test writes into a new run folder under {PATHS['smalltest_dir']} only")
    if not smalltest and real != os.path.realpath(PATHS["out_dir"]):
        stop(f"a full run writes to {PATHS['out_dir']} only")
    os.makedirs(os.path.join(out_dir, "logs"), exist_ok=True)


def check_free_space(path, need_bytes):
    st = os.statvfs(path)
    free = st.f_bavail * st.f_frsize
    if free < MIN_FREE_BYTES + need_bytes:
        stop(f"free space {free:,} bytes on {path} is below {MIN_FREE_BYTES:,} + {need_bytes:,}")
    return free


def md5_of(path, block=1 << 24):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(block), b""):
            h.update(b)
    return h.hexdigest()


def read_sample_metadata():
    """sample_metadata.csv (md5 pinned): 163 rows, one per atlas library, keyed by library_id and read without a row
    filter."""
    got = md5_of(PATHS["sample_metadata"])
    if got != SAMPLE_METADATA_MD5:
        stop(f"sample_metadata.csv md5 {got} differs from the pinned {SAMPLE_METADATA_MD5}")
    sm = pd.read_csv(PATHS["sample_metadata"], dtype=str, keep_default_na=False)
    for c in ("included", "exclusion_reason"):
        if c in sm.columns:
            stop(f"column '{c}' present: this is not the expected sample_metadata.csv layout")
    if len(sm) != N_LIBRARIES or not sm["library_id"].is_unique:
        stop(f"sample_metadata.csv has {len(sm)} rows (expected {N_LIBRARIES} unique libraries)")
    return sm.set_index("library_id", drop=False)


N_SMALLTEST_LIBRARIES = 2  # the S5 small test converts two libraries


def parse_libraries(arg, known):
    """Small test: exactly two comma-separated library ids (the S5 small test converts two
    libraries; they are named with the converter's --libraries option, so there is no default)."""
    if not arg:
        return None
    libs = [x.strip() for x in arg.split(",") if x.strip()]
    bad = [x for x in libs if not re.fullmatch(r"L\d{3}", x) or x not in known]
    if bad or len(set(libs)) != len(libs) or len(libs) != N_SMALLTEST_LIBRARIES:
        stop(f"the small test takes {N_SMALLTEST_LIBRARIES} different known library ids; got {libs} (unknown: {bad})")
    return sorted(libs)


def safe_write_csv(df, path):
    """Write a new CSV (never overwrite): temp file, read-back check, then rename."""
    if os.path.exists(path):
        stop(f"{path} exists; outputs are never overwritten (use a new output folder)")
    tmp = path + ".partial"
    if os.path.exists(tmp):
        stop(f"{tmp} exists from an earlier attempt; it is not removed by this script")
    with open(tmp, "x", newline="", encoding="utf-8") as fh:  # 'x': the OS refuses an existing file
        df.to_csv(fh, index=False)
    back = pd.read_csv(tmp, dtype=str, keep_default_na=False)
    if os.path.getsize(tmp) == 0 or len(back) != len(df) or list(back.columns) != [str(c) for c in df.columns]:
        stop(f"write check failed for {tmp}")
    os.rename(tmp, path)


# ----------------------------------------------------------------------------
# Constant schema fields
# ----------------------------------------------------------------------------
CONSTANT_OBS = {
    "sex_ontology_term_id": ("PATO:0000384", "All donors are male (prostate tissue)."),
    "suspension_type": ("cell", "All libraries are droplet-based 10x Genomics single-cell libraries."),
    "tissue_type": ("tissue", "Primary tissue samples (no cell lines, organoids or cultures)."),
    # False (the author's decision). The note says that every cell comes from the source studies' raw reads,
    # reprocessed with Cell Ranger. In v35d those reads are open (132 libraries), controlled-access (28) or
    # author-provided (3) (access_type). The Luminal file's cells are rows of the all-cells file.
    "is_primary_data": (False, "All cells come from the source studies' raw reads, reprocessed with Cell Ranger, so the "
                               "primary observations belong to the source studies."),
}


def is_primary_data_value():
    """The decided is_primary_data value, which must be a Python bool (a string
    such as 'False' would otherwise be written as True everywhere)."""
    v = CONSTANT_OBS["is_primary_data"]
    if _is_pending(v) or not isinstance(v, tuple) or len(v) != 2 or not isinstance(v[0], bool):
        stop(f"CONSTANT_OBS['is_primary_data'] must be (True or False, note); found {v!r}")
    return v[0]


ORGANISM_ONTOLOGY_TERM_ID = "NCBITaxon:9606"
SCHEMA_OBS = ["assay_ontology_term_id", "cell_type_ontology_term_id", "development_stage_ontology_term_id",
              "disease_ontology_term_id", "donor_id", "is_primary_data", "self_reported_ethnicity_ontology_term_id",
              "sex_ontology_term_id", "suspension_type", "tissue_ontology_term_id", "tissue_type"]
# names CELLxGENE Discover adds itself (as in the earlier verification script); never written by the converter
RESERVED_OBS = ["assay", "cell_type", "development_stage", "disease", "self_reported_ethnicity", "sex", "tissue",
                "observation_joinid", "organism", "organism_ontology_term_id", "ethnicity", "ethnicity_ontology_term_id",
                # schema 7.1.0 (cellxgene-schema 7.1.1 schema_definition.yaml, md5
                # 57a8703d...) reserves experimental_condition and perturbation_types and defines the perturbation
                # fields genetic_perturbation_id, genetic_perturbation_strategy, experimental_condition_ontology_term_id
                "experimental_condition", "perturbation_types", "genetic_perturbation_id", "genetic_perturbation_strategy",
                "experimental_condition_ontology_term_id"]
RESERVED_VAR = ["feature_name", "feature_reference", "feature_biotype", "feature_length", "feature_type"]
# the uns names of schema 7.1.0 (same file): reserved keys, deprecated keys and the schema's own keys. No kept S4 uns
# key may take one of them or the name of the converter's record 'cellxgene_conversion' (check_uns_rule_names)
RESERVED_UNS = ["schema_version", "schema_reference", "citation", "is_pre_analysis"]
DEPRECATED_UNS = ["X_normalization", "default_field", "layer_descriptions", "tags", "version", "contributors",
                  "preprint_doi", "project_description", "project_links", "project_name", "publication_doi"]
SCHEMA_UNS = ["organism_ontology_term_id", "title", "batch_condition", "default_embedding", "X_approximate_distribution",
              "genetic_perturbations"]
CONVERSION_UNS_KEY = "cellxgene_conversion"  # written by convert_to_cellxgene.py

# ----------------------------------------------------------------------------
# S4 object layout, written from the two S4 objects (attributes and small datasets only). Names given by a plan or
# decision source are checked to exist at run time.
# ----------------------------------------------------------------------------
S4_LIBRARY_COLUMN = "library_id"
S4_CELL_TYPE_COLUMN = "cell_type"
# the S4 var columns are gene_ids, feature_types, genome and mt in both objects. gene_ids holds 36,601 unique
# Ensembl gene ids and the var index 36,601 unique gene symbols
S4_GENE_ID_SOURCE = "gene_ids"     # "index" or the S4 var column with the Cell Ranger Ensembl gene ids
S4_GENE_SYMBOL_SOURCE = "index"    # "index" or the S4 var column with the gene symbols
# The rule sets below (obs, var, layers, uns) are complete for the S4 objects. Each slot's constant is
# "<record>; <md5>". The md5 is that of the listing (object, file, bytes, md5) of the two S4 objects
# that the rules were written from. The 'keep' rules of 7 obs columns and of the uns key 'embeddings' need
# the author's confirmation (OBS_KEEP_CONFIRMED below).
S4_OBS_RULES_FROM = "S4 object listing; 7c15f97ac0a28a1710a834d43e9a0455"
S4_VAR_RULES_FROM = "S4 object listing; 7c15f97ac0a28a1710a834d43e9a0455"
S4_LAYER_RULES_FROM = "S4 object listing; 7c15f97ac0a28a1710a834d43e9a0455"
S4_UNS_RULES_FROM = "S4 object listing; 7c15f97ac0a28a1710a834d43e9a0455"
OBS_ACTIONS = ("keep", "rename", "drop")
VAR_ACTIONS = ("keep", "rename", "drop")
LAYER_ACTIONS = ("drop",)        # layers are not written
UNS_ACTIONS = ("keep", "drop")   # any other word stops the run

# obs rule: (action, new name, note for cellxgene_field_mapping.csv); action keep | rename | drop.
# An S4 obs column named like a sample_metadata.csv column needs no entry: it is kept and checked equal to
# sample_metadata.csv per library as exact text, also with an entry here. To write
# the sample_metadata.csv value instead, give it a 'drop' entry: that field is then added (SAMPLE_FIELDS_ADD), and
# cellxgene_field_mapping.csv lists it as added. A rename into a sample_metadata.csv
# column name stops. Every other S4 obs column needs an entry, and a column without one
# stops the run. The rules cover every obs column of the S4 objects.
S4_OBS_RULES = {
    # sample_metadata.csv has a column of this name. So the converter also checks it equal to
    # sample_metadata.csv per library as exact text, and it adds no second library column.
    S4_LIBRARY_COLUMN: ("keep", S4_LIBRARY_COLUMN,
                        "Library identifier (Lnnn), the key to join sample_metadata.csv."),
    # reserved schema names avoided; author annotations retained (article, Table 2)
    S4_CELL_TYPE_COLUMN: ("rename", "author_cell_type",
                  "Atlas cell type, renamed because the schema reserves 'cell_type' for the Cell Ontology label that "
                  "CELLxGENE Discover adds. Source of cell_type_ontology_term_id."),
    # its labels as the S4 objects hold them: the 18 atlas cell types, with Dendritic_plasmacytoid (within
    # Dendritic) and Lymphatic_endothelial (within Endothelial cell)
    "cell_type_fine": ("rename", "author_cell_type_fine",
                       "Finer atlas cell type, renamed to match author_cell_type. Its labels are the 18 atlas cell types, "
                       "with Dendritic_plasmacytoid and Lymphatic_endothelial split from Dendritic and Endothelial cell. "
                       "In luminal_cells_cellxgene_v1.h5ad, its only label is Luminal."),
    # column names luminal_subgroup and luminal_label_flag
    # notes without a metadata_dictionary.csv pointer: the v35d dictionary does not cover these obs columns
    "luminal_subgroup": ("keep", "luminal_subgroup",
                         "Putative luminal subgroup label (Methods), missing for cells that are not Luminal."),
    "luminal_label_flag": ("keep", "luminal_label_flag",
                           "Reason to use the luminal subgroup label with caution (Usage Notes), missing where there is "
                           "no flag."),
    # the 7 other S4 obs columns, each 'keep', so that the CELLxGENE files hold the fields of the S4 objects. A full run
    # stops until the author confirms them (OBS_KEEP_CONFIRMED). A small test does not check this. dtypes as in the S4
    # objects (int32, int32, float64, float64, bool, category, category). keep_cat writes the '' of
    # first_submission_cell_id (cells that are not first-submission cells) as a missing value, which same_column accepts.
    "nCount_RNA_allgenes": ("keep", "nCount_RNA_allgenes",
                            "Number of UMIs of the cell over all 36,601 genes (13 mitochondrial genes included)."),
    "nFeature_RNA_allgenes": ("keep", "nFeature_RNA_allgenes",
                              "Number of genes with at least one UMI in the cell, over all 36,601 genes."),
    "pct_mt": ("keep", "pct_mt", "Percentage of the cell's UMIs from the 13 mitochondrial genes."),
    "nuclear_fraction_pct": ("keep", "nuclear_fraction_pct",
                             "Intronic UMI fraction of the cell, in percent."),
    "luminal_object": ("keep", "luminal_object",
                       "True for the cells of the luminal integration (the cells of luminal_cells_cellxgene_v1.h5ad), "
                       "False for every other cell."),
    "luminal_training_class": ("keep", "luminal_training_class",
                               "Luminal training class of the cell in the luminal scANVI training (Benign, BPH, ERG, ETV, "
                               "nonETS, AR or NEPC). Luminal-integration cells without a class are 'unlabeled', and cells "
                               "outside that integration are missing."),
    "first_submission_cell_id": ("keep", "first_submission_cell_id",
                                 # 2,748 (all cells) and 799 (luminal) L138/L139 cells have an obs index ending
                                 # in '-2' and this value ending in '-1'
                                 "Identifier of this cell in the data of the article's first submitted version, for "
                                 "joining those data. Missing for cells that were not in those data. Where "
                                 "present, it equals the obs index, except that it ends in '-1' where the obs index ends "
                                 "in '-2' (L138 and L139)."),
}
# var rule: (action, new name, note). The gene id and symbol sources are handled separately.
S4_VAR_RULES = {
    # the 13 mitochondrial genes stay in the released raw counts with var['mt']
    "mt": ("keep", "mt", "Mitochondrial gene flag (13 genes). These genes are kept in the counts."),
    # Cell Ranger constant columns, as in the earlier conversion; the constant is checked before the drop
    "feature_types": ("drop", "", "Constant 'Gene Expression'. CELLxGENE Discover adds feature_type and feature_biotype itself."),
    "genome": ("drop", "", "Constant 'GRCh38'. CELLxGENE Discover adds feature_reference itself."),
    # no other S4 var column (gene_ids is S4_GENE_ID_SOURCE, the index S4_GENE_SYMBOL_SOURCE)
}
S4_VAR_DROP_IF_CONSTANT = {"feature_types": "Gene Expression", "genome": "GRCh38"}
# gene_id_remapping.csv column highly_variable_all_cells (the author's decision): "take the HVG column from S4 if S4
# has one all-cells HVG flag; otherwise omit it" (the earlier conversion took it from the first-version var['highly_variable']).
# S4_ALL_CELLS_HVG_FLAGS is the S4 fact this needs: the var columns of the S4 all-cells object that hold an
# all-cells HVG flag, written from the S4 var listing together with S4_VAR_RULES (S4_VAR_RULES_FROM); [] if there
# is none. remap_hvg_column() applies that decision to it.
# The S4 var has no HVG flag column (its only boolean column is mt).
S4_ALL_CELLS_HVG_FLAGS = []
REMAP_HVG_OUT_COLUMN = "highly_variable_all_cells"  # the gene_id_remapping.csv column name (as before)


def remap_hvg_column():
    """The S4 all-cells HVG flag column when S4 has exactly one, else None (the column
    highly_variable_all_cells is left out of gene_id_remapping.csv)."""
    flags = S4_ALL_CELLS_HVG_FLAGS
    if _is_pending(flags) or not DECIDED_FORMS["S4_ALL_CELLS_HVG_FLAGS"][0](flags):
        stop(f"S4_ALL_CELLS_HVG_FLAGS must be {DECIDED_FORMS['S4_ALL_CELLS_HVG_FLAGS'][1]}; found {flags!r}")
    return flags[0] if len(flags) == 1 else None


def hvg_flag_values(var):
    """Values of highly_variable_all_cells from the S4 all-cells var, or None when the column is left out.
    Every column listed in S4_ALL_CELLS_HVG_FLAGS must be in var and be boolean (a flag), else the run stops;
    called before the X pass of build_gene_id_remap.py (not included)."""
    col = remap_hvg_column()
    for c in S4_ALL_CELLS_HVG_FLAGS:
        if c not in var.columns:
            stop(f"S4_ALL_CELLS_HVG_FLAGS names '{c}', which is not a var column of the all-cell source object")
        if not pd.api.types.is_bool_dtype(var[c].dtype):
            stop(f"S4_ALL_CELLS_HVG_FLAGS names '{c}', whose dtype {var[c].dtype} is not boolean (not a flag)")
    return None if col is None else var[col].values


def s4_gene_ids_symbols(var):
    """Cell Ranger Ensembl gene ids and gene symbols of an S4 var table (sources set above)."""
    def get(src, what):
        if src == "index":
            return var.index.astype(str).values
        if src not in var.columns:
            stop(f"var column '{src}' ({what}) is not in the source object")
        return var[src].astype(str).values
    ids, sym = get(S4_GENE_ID_SOURCE, "gene id"), get(S4_GENE_SYMBOL_SOURCE, "gene symbol")
    if len(ids) != N_RELEASE_GENES or not pd.Index(ids).is_unique:
        stop(f"{len(ids)} gene ids (expected {N_RELEASE_GENES} unique)")
    return ids, sym


# the S4 objects have no layer (X holds the raw counts), so the layer rules are empty and complete. Their uns has one
# key, 'embeddings', with one entry per obsm key. Each entry gives method, seed, coordinates, dimensions, cell labels
# and software, plus a note or the training labels where the S4 script wrote one.
S4_LAYER_RULES = {}  # one (action, note) per S4 layer
S4_UNS_RULES = {     # one (action, note) per S4 uns key that is not '<column>_colors'
    # 'keep' after the author's confirmation (OBS_KEEP_CONFIRMED): every obsm embedding is copied, and this record
    # describes each of them
    "embeddings": ("keep", "Description of each obsm embedding: method, seed, coordinates, dimensions, cell_labels "
                           "(label-free or label-guided) and software. scANVI entries add training_labels (atlas cell "
                           "types or luminal training classes), and BBKNN and scANVI entries add a note. Copied unchanged, except "
                           "that the all-cell training_labels name author_cell_type."),
}
# The training_labels of uns['embeddings'] name, in parentheses, the S4 obs column of the training labels.
# In the CELLxGENE files, the all-cell atlas labels are written as author_cell_type (S4_OBS_RULES).
# So the converter rewrites that one token and nothing else. luminal_training_class keeps its name.
# Every training_labels text must hold its file's token exactly once, and each file must have some.
TRAINING_LABEL_COLUMN = {"all_cells": S4_CELL_TYPE_COLUMN, "luminal": "luminal_training_class"}


def rename_training_label_tokens(embeddings, dataset):
    """uns['embeddings'] with '(<S4 column>)' in each training_labels rewritten to '(<written column>)'; stops if a
    training_labels text does not hold that token exactly once, or if no entry has training_labels."""
    col = TRAINING_LABEL_COLUMN[dataset]
    old, new = f"({col})", f"({S4_OBS_RULES[col][1]})"
    out, n = {}, 0
    if not isinstance(embeddings, dict):  # a stop instead of an AttributeError
        stop(f"uns['embeddings'] of the {dataset} source object is a {type(embeddings).__name__}, not a dict")
    for key, entry in embeddings.items():
        if not isinstance(entry, dict):  # a stop instead of a TypeError
            stop(f"uns['embeddings']['{key}'] of the {dataset} source object is a {type(entry).__name__}, not a dict")
        entry = dict(entry)
        if "training_labels" in entry:
            text = str(entry["training_labels"])
            if text.count(old) != 1:
                stop(f"uns['embeddings']['{key}'] training_labels should name {old} exactly once: {text!r}")
            entry["training_labels"] = text.replace(old, new)
            n += 1
        out[key] = entry
    if n == 0:
        stop(f"uns['embeddings'] of the {dataset} source object has no training_labels (expected the scANVI entries)")
    return out


# The author's confirmation of the 'keep' rules above (7 S4 obs columns of S4_OBS_RULES and the uns key 'embeddings').
# It is PENDING until the author confirms, then a record of the confirmation. A full run stops while it is PENDING
# (check_full_run_confirmed). A small test does not check it.
OBS_KEEP_CONFIRMED = "confirmed by the author"


def check_uns_rule_names():
    """Stop when S4_UNS_RULES keeps an S4 uns key under a reserved, deprecated or schema uns name of schema 7.1.0, or
    under the converter's record name. Reads the rules only."""
    taken = set(RESERVED_UNS) | set(DEPRECATED_UNS) | set(SCHEMA_UNS) | {CONVERSION_UNS_KEY}
    bad = sorted(k for k, r in S4_UNS_RULES.items() if r[0] == "keep" and k in taken)
    if bad:
        stop(f"source uns keys kept under a reserved, deprecated or schema uns name (schema {SCHEMA_VERSION}): {bad}")


def check_full_run_confirmed():
    """Full run only: stop while OBS_KEEP_CONFIRMED is PENDING or of the wrong form. Makes no folder or file."""
    require_decided(FULL_RUN_NEEDED)
    check_forms(FULL_RUN_NEEDED)
    print(f"full-run confirmation present: OBS_KEEP_CONFIRMED = {OBS_KEEP_CONFIRMED}", flush=True)

# ----------------------------------------------------------------------------
# Sample-level fields added to the CELLxGENE obs from sample_metadata.csv (values per library)
# ----------------------------------------------------------------------------
# the earlier sample-level obs columns in the v33 names (prostate_weight_g no longer exists). Not used while
# SAMPLE_FIELDS_ADD is 'all'. v35d lacks three of these names (internal_sample_id, sample_label and study_short_id),
# so this list would stop as missing fields
SAMPLE_FIELDS_V2_SET = [
    "internal_sample_id", "sample_id", "sample_label", "patient_id", "study_id", "study_short_id", "study_patient_id",
    "run_accessions", "disease_grade_category", "gleason_score", "gleason_score_at_diagnosis", "biopsy_grade_group",
    "tissue_site", "tissue_site_group", "surgical_procedure", "surgical_margin", "tnm_stage", "tnm_stage_at_mcrpc",
    "age", "race", "bmi", "hypertension", "hyperlipidemia", "diabetes", "mri_findings", "psa_at_diagnosis_ng_ml",
    "psa_at_sampling_ng_ml", "pretreatment_adt", "first_line_therapy", "second_line_therapy",
    "AMACR_reported", "TP63_reported", "AR_reported", "TP53_reported", "RB1_reported", "ERG_reported", "PTEN_reported",
    "AKT1_reported", "CDKN2A_reported", "CDKN2B_reported", "MYC_reported",
    "sample_source_institute", "sample_source_city", "sample_source_country", "histology_image_in_source", "curation_notes",
    "disease_state", "disease_subtype", "gleason_grade_group", "library_chemistry", "facs_sorted_study", "sorted_fraction",
    "study_accession", "sample_accession",
]
# "all" (the author's decision): every column of sample_metadata.csv
# (the article, Table 2: curated sample metadata as extra columns); 60 columns in v35d.
SAMPLE_FIELDS_ADD = "all"
N_SAMPLE_COLUMNS = 60  # v35d header; 'all' stops on another count
NUMERIC_SAMPLE_FIELDS = {"bmi"}  # as before (stored as float); every other field is categorical text
REQUIRED_SAMPLE_FIELDS = ["study_id"]  # uns['batch_condition'] names it


def sample_fields_to_add(sm_columns):
    if SAMPLE_FIELDS_ADD == "all":
        fields = list(sm_columns)
        if len(fields) != N_SAMPLE_COLUMNS or len(set(fields)) != len(fields):
            stop(f"SAMPLE_FIELDS_ADD = 'all' expects the {N_SAMPLE_COLUMNS} distinct sample_metadata.csv columns; "
                 f"found {len(fields)} ({len(set(fields))} distinct)")
    else:
        fields = list(SAMPLE_FIELDS_ADD)
    missing = [c for c in fields + REQUIRED_SAMPLE_FIELDS if c not in sm_columns]
    if missing:
        stop(f"sample fields not in sample_metadata.csv: {missing}")
    for c in REQUIRED_SAMPLE_FIELDS:
        if c not in fields:
            fields.append(c)
    return fields


def obs_rule(col, sample_columns):
    """(action, new name, note) for an S4 obs column, or None when its rule is PENDING."""
    if col in S4_OBS_RULES:
        return S4_OBS_RULES[col]
    if col in sample_columns:
        return ("keep", col, "Sample-level field (see metadata_dictionary.csv).")
    return None


# ----------------------------------------------------------------------------
# Cell type: the harmonized labels of the S4 obs['cell_type'] (author labels are kept as
# author_cell_type; the CL term is added in cell_type_ontology_term_id)
# 'Mast cell' CL:0000097 and 'Neutrophil' CL:0000775 added,
# 'Granulocytes' dropped, the Monocyte note corrected (neutrophils are now their own label).
# ----------------------------------------------------------------------------
# columns: author label -> (CL id, match type, rationale, needs author review?)
CELL_TYPE_MAP = {
    "Luminal": ("CL:0002340", "exact",
                "Luminal epithelial cells of the prostate. Also used for luminal (tumor) cells found in metastatic sites, which are of prostate luminal origin; malignant status is not encoded in the cell type.", "yes"),
    "Basal": ("CL:0002341", "exact", "Basal epithelial cells of the prostate (all cells come from prostate tissue).", "no"),
    "Club": ("CL:4032000", "close",
             "CL term 'club-like cell of the urethral epithelium' (synonym 'prostate club cell'; SCGB1A1+), defined from the human prostate/prostatic urethra (Henry et al. 2018). Preferred over CL:0000158 'club cell', which is a tracheobronchial (lung) cell type.", "yes"),
    "Hillock": ("CL:0017011", "close",
                "CL term 'hillock cell of prostatic urethral epithelium' (KRT13+), defined from the human prostate (Henry et al. 2018). Generic alternative: CL:4030024 'hillock cell'.", "yes"),
    "Ionocyte": ("CL:0005006", "broader",
                 "Generic 'ionocyte'; CL has no prostate-specific ionocyte term (only pulmonary and salivary gland ionocytes).", "no"),
    "T cell": ("CL:0000084", "exact", "", "no"),
    "NK cell": ("CL:0000623", "exact", "natural killer cell", "no"),
    "B cell": ("CL:0000236", "exact", "", "no"),
    "Plasma": ("CL:0000786", "exact", "plasma cell", "no"),
    "Macrophage": ("CL:0000235", "exact", "", "no"),
    # review flag back to the earlier 'yes' (no row decides 'no'). The earlier note described the
    # first-version objects (498 cells labeled 'Neutrophils' at the finer level).
    # MONOCYTE_S4_COUNT_FROM records how many of them keep the S4 'Monocyte' label.
    # The flag must stay 'yes' while that count is PENDING. check_review_flags() enforces this, and
    # check_handover() blocks the hand-over until the count is recorded.
    "Monocyte": ("CL:0000576", "exact", "monocyte", "yes"),
    "Dendritic": ("CL:0000451", "exact",
                  "dendritic cell (group contains conventional and plasmacytoid DCs at the finer level; parent term used).", "no"),
    "Mast cell": ("CL:0000097", "exact", "mast cell", "no"),
    "Neutrophil": ("CL:0000775", "exact", "neutrophil", "no"),
    "Endothelial cell": ("CL:0000115", "exact", "endothelial cell (group also contains lymphatic EC and a small adipocyte-like subset at the finer level).", "no"),
    "Fibroblast": ("CL:0000057", "exact",
                   "Generic term used for all tissues (cells come from prostate and metastatic sites). Prostate-specific alternative for prostate tissue: CL:1000299 'fibroblast of connective tissue of prostate'.", "no"),
    "SMC": ("CL:0000192", "exact",
            "smooth muscle cell. Prostate-specific alternative for prostate tissue: CL:1000487 'smooth muscle cell of prostate' (not used because the cluster also contains vascular smooth muscle/pericyte-like cells and cells from metastatic sites).", "no"),
    "Schwann": ("CL:0002573", "exact", "Schwann cell", "no"),
}

# ----------------------------------------------------------------------------
# Assay: sample_metadata.csv 'library_chemistry'; L083 by library (replaces the
# earlier EFO:0009922)
# ----------------------------------------------------------------------------
ASSAY_MAP = {
    "10x 3' v2": ("EFO:0009899", "Schema-recommended term for 10x 3' v2."),
    "10x 3' v3": ("EFO:0009922", "Schema-recommended term for 10x 3' v3."),
    "10x 3' v3.1": ("EFO:0022980", "EFO term '10x 3' v3.1' (child of EFO:0009922 '10x 3' v3')."),
    "10x 5' v1": ("EFO:0011025", "EFO term '10x 5' v1'."),
    "10x 5' v2": ("EFO:0009900", "Schema-recommended term for 10x 5' v2."),
    "10x 3' (version not reported)": ("EFO:0030003",
                                       "Chemistry version not reported by the source study (P21, GSE206962): parent term '10x 3' transcription profiling'. Can be refined if the Cell Ranger auto-detected chemistry is checked."),
    # No entry for an empty library_chemistry (the author's decision). The v35d dictionary allows it
    # ("empty = the source does not state the version and the Cell Ranger detection cannot decide it, or the source
    # conflicts with the data"); in v35d only L083 is empty, and L083 has
    # its decided library term (ASSAY_BY_LIBRARY), which assay_term() takes first. An empty value
    # of any other library stops as an unknown key and comes back to the author. Never filled by inference.
}
ASSAY_BY_LIBRARY = {
    "L083": ("EFO:0030003",
             "The source gives the 10x 3' v2 read layout, while Cell Ranger detects 'Single Cell 3' v3' (which it cannot tell from v3.1): parent term '10x 3' transcription profiling'."),
}

# ----------------------------------------------------------------------------
# Disease. Non-metastatic groups by disease_state, as in the earlier conversions.
# The mCRPC groups by the two status columns of sample_metadata.csv, keyed on the full
# strings, plus a subtype term for the 8 named NEPC libraries
# (MONDO:0006390 for the 6 small-cell libraries; MONDO:0002477 for L079 and L082,
# which replaces an earlier term). Multi-term values are sorted and de-duplicated
# (schema 7.1.0 multi_term sorted). Nothing falls back: an unknown key stops the run.
# ----------------------------------------------------------------------------
DISEASE_BY_STATE = {
    "normal": ("PATO:0000461", "Organ donors / non-cancer prostates: 'normal'.", "no"),
    "benign": ("MONDO:0005082",
               "Histologically benign tissue adjacent to tumor, collected from prostate adenocarcinoma patients (radical prostatectomy). Annotated with the donor's disease; the sample-level state is kept in obs['disease_state'] = 'benign'. Alternative: 'PATO:0000461' (normal) if the authors prefer a tissue-level convention.", "yes"),
    "BPH": ("MONDO:0010811", "benign prostatic hyperplasia", "no"),
    "primary_tumor": ("MONDO:0005082",
                      "Primary prostate adenocarcinoma; the libraries with a disease_subtype value take the same term.", "no"),
}
MCRPC_STATES = ("mCRPC", "mCRPC-NEPC")  # the two mCRPC group names
# author_review_needed of the mCRPC rows of disease_ontology_mapping.csv: 'no', because the author decided these
# terms (terms from the status columns; MONDO:0006390 for the 6 small-cell NEPC libraries; and
# MONDO:0002477 for L079 and L082). The earlier conversion: 'no' for mCRPC, 'yes' for its teSCC and NETD rows.
MCRPC_AUTHOR_REVIEW = "no"
DISEASE_BY_STATUS = {
    # (castration_status, metastatic_status) -> (terms, rationale)
    ("castration-resistant", "yes"): (("MONDO:0004956", "MONDO:0850353"),
                                      "Castration-resistant and metastatic: 'metastatic prostate carcinoma' and 'castration-resistant prostate carcinoma'."),
    ("castration-resistant", "no"): (("MONDO:0850353",),
                                     "Castration-resistant, not metastatic per the status columns: 'castration-resistant prostate carcinoma'."),
    ("castration-resistant", "not reported"): (("MONDO:0850353",),
                                               "Castration-resistant, metastatic status not reported: 'castration-resistant prostate carcinoma'."),
    ("not determinable from the source", "yes"): (("MONDO:0004956",),
                                                  "Metastatic; castration status not determinable from the source: 'metastatic prostate carcinoma'."),
    ("conflicting in the source (called CRPC, no ADT listed)", "yes"): (("MONDO:0004956",),
                                                                        "Metastatic; castration status conflicting in the source: 'metastatic prostate carcinoma'."),
}
# the decided counts (mCRPC + mCRPC-NEPC), checked at run time
EXPECTED_STATUS_COUNTS = {
    ("castration-resistant", "yes"): 20,
    ("castration-resistant", "no"): 2,
    ("castration-resistant", "not reported"): 2,
    ("not determinable from the source", "yes"): 9,
    ("conflicting in the source (called CRPC, no ADT listed)", "yes"): 1,
}
N_MCRPC_LIBRARIES = 34  # 26 mCRPC + 8 mCRPC-NEPC
NEPC_SUBTYPE_TERM = {
    **{lib: ("MONDO:0006390", "Small-cell neuroendocrine subtype: adds 'prostate small cell carcinoma'.")
       for lib in ("L064", "L072", "L073", "L083", "L085", "L086")},
    **{lib: ("MONDO:0002477", "Well-differentiated neuroendocrine subtype: adds 'prostate neuroendocrine neoplasm'.")
       for lib in ("L079", "L082")},
}
EXPECTED_DISEASE = {
    "L079": "MONDO:0002477 || MONDO:0850353",
    "L082": "MONDO:0002477 || MONDO:0004956 || MONDO:0850353",
}
# Cross-checks of the by-library lists. The strings were quoted
# from the v33d dictionary and re-checked against the v35d dictionary and values.
# A mismatch stops the run and never changes a term.
NEPC_EXPECTED_NEUROENDOCRINE_TYPE = {"MONDO:0006390": "small cell", "MONDO:0002477": "well-differentiated"}
WONG_STUDY_ID = "P02_2022_NatCommun_Wong"  # the study of the 14 libraries of DEV_STAGE_BY_LIBRARY
N_PRIMARY_SUBTYPE_LIBRARIES = 11  # the 11 primary subtype libraries keep MONDO:0005082

# ----------------------------------------------------------------------------
# Tissue: sample_metadata.csv 'tissue_site' (+ zone in curation_notes for prostate), with the
# decided library terms first
# ----------------------------------------------------------------------------
TISSUE_SITE_MAP = {
    "Prostate": ("UBERON:0002367", "prostate gland (zone not reported)", "no"),
    "Bone": ("UBERON:0002371",
             "Kfoury et al. (GSE143791) GEO records: 'tissue: bone marrow', 'tissue type: solid metastatic tissue (Tumor fraction)' from vertebral bone metastases (spinal cord compression surgery). Alternative: UBERON:0002481 'bone tissue'.", "yes"),
    "Liver": ("UBERON:0002107", "liver", "no"),
    "Brain": ("UBERON:0000955", "brain", "no"),
    "Epidural": ("UBERON:0003691", "epidural space (epidural metastasis)", "no"),
    "Mediastinal mass": ("UBERON:0003728", "mediastinum (mediastinal mass)", "no"),
    "SupraclavicularLymphNode": ("UBERON:0035279", "supraclavicular lymph node", "no"),
    "PosteriormedialExternalLymphNode": ("UBERON:0015880",
                                         "external iliac lymph node ('posteromedial external' node in the source table; UBERON has medial/intermediate/lateral external iliac nodes only, so the parent term is used).", "yes"),
    "RetrocruralAbdomenMass": ("UBERON:0000916",
                               "abdomen (retrocrural abdominal mass; not a lymph node per the source label). Alternative if the mass was nodal: UBERON:0038054 'retrocrural lymph node'.", "yes"),
    "ParaAorticLymphNode": ("UBERON:0016386", "paraaortic lymph node", "no"),
    "PelvicLymphNode": ("UBERON:0015876", "pelvic lymph node", "no"),
    "PelvicSideWall": ("UBERON:8480043", "pelvic wall (pelvic side-wall metastasis)", "no"),
    # No site rule for the two compound values (the author's decision): in v35d,
    # 'brain or epidural (library not linked to biopsy)' is held by L067, L068 and L071 only, and
    # 'spinal epidural tumor (paper); bone marrow (GEO)' by L101-L109 only; all twelve have
    # their decided library term (TISSUE_BY_LIBRARY), which tissue_term() takes first.
    # Either value held by any other library stops as an unknown key and comes back to the author. Never filled by inference.
}
# author_review_needed of tissue rows whose term is set by library (other than prostate zone rows): 'no',
# because the author decided these terms (UBERON:0010000 for the mixed-tissue libraries, source text
# kept, and L101-L109 UBERON:0001434). The earlier conversion had no rule by library.
TISSUE_LIBRARY_AUTHOR_REVIEW = "no"
# The libraries those two rows cover: the author's decision counts 12 mixed-tissue libraries (the count, not the ids; the
# ids are those of the TISSUE_BY_LIBRARY entries below), of which a later decision sets L101-L109, leaving L067,
# L068 and L071. Only a row made of these libraries takes
# TISSUE_LIBRARY_AUTHOR_REVIEW. L146 and L151 are set by library, not by these rows; tissue_term
# stops for them at a site other than 'Prostate'.
TISSUE_LIBRARY_DECIDED = frozenset(["L067", "L068", "L071"] + [f"L{i:03d}" for i in range(101, 110)])
TISSUE_BY_LIBRARY = {
    **{lib: ("UBERON:0010000",
             "The source names two sites for this library; 'multicellular anatomical structure' is used, and the source text is kept in tissue_site.")
       for lib in ("L067", "L068", "L071")},
    **{f"L{i:03d}": ("UBERON:0001434", "'skeletal system'; the source text is kept in tissue_site.")
       for i in range(101, 110)},
    **{lib: ("UBERON:8410026", "peripheral zone of prostate (zone reported by the source study)")
       for lib in ("L146", "L151")},  # by name (the earlier zone regex misses their notes)
}
PROSTATE_ZONE_RULES = [
    # (regex on curation_notes, term, rationale); peripheral rule widened (earlier: ^PeripheralZone|^Peripheral$).
    # v35d writes the first curation_notes entry of 20 Prostate libraries as 'zone: <zone> (<source>)' and the bare
    # '<zone>' for 8 others. So each rule allows an optional 'zone: ' label before the zone.
    (r"^(?:zone: )?PeripheralZone|^(?:zone: )?Peripheral\b", "UBERON:8410026", "peripheral zone of prostate (zone reported by the source study; curation_notes)"),
    (r"^(?:zone: )?TransitionZone", "UBERON:8410025", "transition zone of prostate (zone reported by the source study; curation_notes; Henry et al. transition-zone samples include the prostatic urethra)"),
]

# ----------------------------------------------------------------------------
# Self-reported ethnicity: sample_metadata.csv 'race' (all reporting cohorts are US-based)
# ----------------------------------------------------------------------------
ETHNICITY_MAP = {
    "White": ("HANCESTRO:0590", "European American (HANCESTRO synonym 'White'; US cohorts)."),
    "Black or African American": ("HANCESTRO:0568", "African American (US cohort)."),
    "Asian": ("HANCESTRO:0847", "Asian (geography-based population category)."),
    "Hispanic": ("HANCESTRO:0612", "Hispanic or Latin."),
    "6 White and 1 Black (African American)": ("unknown", "Cohort-level statement only (Wong et al.); individual donors cannot be assigned."),
    "": ("unknown", "Not reported."),
    # The v35d dictionary allows 'Asian | Black or African American | Hispanic | White' and empty ("empty where no race
    # or ethnicity is stated for the donor"). The v35d race values are '', 'Asian', 'Black or African American',
    # 'Hispanic' and 'White'; each has its earlier rule above.
    # The cohort-level key '6 White and 1 Black (African American)' is held by no v35d library (an earlier rule, kept).
    # Any other value stops as an unknown key and comes back to the author.
}

# ----------------------------------------------------------------------------
# Development stage from sample_metadata.csv 'age'; the 14 Wong libraries by library
# ----------------------------------------------------------------------------
DEV_STAGE_BY_LIBRARY = {
    f"L{i:03d}": ("HsapDv:0000258",
                  "Cohort-level age only (Wong et al., mean 65, range 52-74): all donors are adults, so 'adult stage' is used.")
    for i in range(1, 15)
}


def age_to_hsapdv(age, label_to_id):
    """Return (term_id, rule) for a verbatim age string; an unreadable age stops the run."""
    a = str(age).strip()
    # the earlier not-reported strings; any other string that is not a number stops. In v35d the
    # only age that is not a number is '' (dictionary: "empty where no age is stated for the donor").
    if a in ("", "NA", "nan"):
        return "unknown", "Age not reported by the source study."
    try:
        years = float(a)
    except ValueError:
        stop(f"unreadable age '{a}' (no rule)")
    if not np.isfinite(years):  # 'NaN', 'inf': a STOP, not an error from int()
        stop(f"unreadable age '{a}' (no rule)")
    if years != int(years):
        stop(f"non-integer age {a}")
    years = int(years)
    tid = label_to_id.get(f"{years}-year-old stage")
    if tid is None:
        stop(f"no HsapDv term for age {years}")
    return tid, f"Age in years reported per donor -> HsapDv '{years}-year-old stage'."


def sorted_terms(terms):
    """Schema 7.1.0 multi-term value: unique terms, sorted, joined by ' || '."""
    return " || ".join(sorted(set(terms)))


def assay_term(r):
    lib = r["library_id"]
    if lib in ASSAY_BY_LIBRARY:
        return ASSAY_BY_LIBRARY[lib]
    chem = r["library_chemistry"]
    if chem not in ASSAY_MAP or _is_pending(ASSAY_MAP[chem]):
        stop(f"{lib}: library_chemistry '{chem}' has no decided assay rule")
    return ASSAY_MAP[chem]


def disease_term(r):
    lib, state = r["library_id"], r["disease_state"]
    if state in DISEASE_BY_STATE:
        if lib in NEPC_SUBTYPE_TERM:
            stop(f"{lib}: NEPC subtype library with disease_state '{state}'")
        tid, why, _ = DISEASE_BY_STATE[state]
        return tid, why
    if state not in MCRPC_STATES:
        stop(f"{lib}: disease_state '{state}' has no disease rule")
    key = (r["castration_status"], r["metastatic_status"])
    if key not in DISEASE_BY_STATUS:
        stop(f"{lib}: status pair {key} has no disease rule")
    terms, why = DISEASE_BY_STATUS[key]
    terms = list(terms)
    if lib in NEPC_SUBTYPE_TERM:
        t, w = NEPC_SUBTYPE_TERM[lib]
        terms.append(t)
        why = why + " " + w
    return sorted_terms(terms), why


def tissue_term(r):
    """Per-library tissue term: decided library term, else tissue_site with the prostate zone rule."""
    lib, site = r["library_id"], r["tissue_site"]
    if lib in ("L146", "L151") and site != "Prostate":  # their library term is a prostate zone
        stop(f"{lib}: tissue_site '{site}' is not 'Prostate', but its library term {TISSUE_BY_LIBRARY[lib][0]} is a prostate zone")
    if lib in TISSUE_BY_LIBRARY:
        return TISSUE_BY_LIBRARY[lib]
    if site == "Prostate":
        notes = str(r.get("curation_notes", ""))
        for rx, tid, why in PROSTATE_ZONE_RULES:
            if re.search(rx, notes):
                return tid, why
    if site not in TISSUE_SITE_MAP or _is_pending(TISSUE_SITE_MAP[site]):
        stop(f"{lib}: tissue_site '{site}' has no decided tissue rule")
    tid, why, _ = TISSUE_SITE_MAP[site]
    return tid, why


def ethnicity_term(r):
    race = r["race"]
    if race not in ETHNICITY_MAP or _is_pending(ETHNICITY_MAP[race]):
        stop(f"{r['library_id']}: race '{race}' has no decided ethnicity rule")
    return ETHNICITY_MAP[race]


def development_stage_term(r, label_to_id):
    lib = r["library_id"]
    if lib in DEV_STAGE_BY_LIBRARY:
        return DEV_STAGE_BY_LIBRARY[lib]
    return age_to_hsapdv(r["age"], label_to_id)


def per_sample_schema_fields(sm, hsapdv_label_to_id):
    """Return one row per library (library_id) with all sample-level schema fields and the
    rule used for each; checks the decided terms (check_decided_terms)."""
    sm = sm.reset_index(drop=True)
    # donor_id = study_patient_id: one study per donor id, else two donors would merge
    n_studies = sm.groupby("study_patient_id")["study_id"].nunique()
    if (n_studies > 1).any():
        stop(f"study_patient_id values used by more than one study_id: {sorted(n_studies[n_studies > 1].index)}")
    rows = []
    for _, r in sm.iterrows():
        assay, assay_rule = assay_term(r)
        dis, dis_rule = disease_term(r)
        dev, dev_rule = development_stage_term(r, hsapdv_label_to_id)
        tis, tis_rule = tissue_term(r)
        eth, eth_rule = ethnicity_term(r)
        rows.append(dict(
            library_id=r["library_id"],
            assay_ontology_term_id=assay,
            disease_ontology_term_id=dis,
            development_stage_ontology_term_id=dev,
            self_reported_ethnicity_ontology_term_id=eth,
            tissue_ontology_term_id=tis,
            donor_id=r["study_patient_id"],
            assay_rule=assay_rule,
            disease_rule=dis_rule,
            development_stage_rule=dev_rule,
            self_reported_ethnicity_rule=eth_rule,
            tissue_rule=tis_rule,
        ))
    per = pd.DataFrame(rows).set_index("library_id")
    check_decided_terms(sm.set_index("library_id", drop=False), per)
    return per


def check_decided_terms(sm, per):
    """Every decided library term and count; any miss stops."""
    errs = []
    for lib in [f"L{i:03d}" for i in range(1, 11)]:
        if lib in per.index and per.loc[lib, "assay_ontology_term_id"] != "EFO:0011025":
            errs.append(f"{lib} assay {per.loc[lib, 'assay_ontology_term_id']} (decided EFO:0011025)")
    fam = sm[sm["disease_state"].isin(MCRPC_STATES)]
    if len(fam) != N_MCRPC_LIBRARIES:
        errs.append(f"{len(fam)} libraries in {MCRPC_STATES} (decided table: {N_MCRPC_LIBRARIES})")
    got = fam.groupby(["castration_status", "metastatic_status"]).size().to_dict()
    if got != EXPECTED_STATUS_COUNTS:
        errs.append(f"status-pair counts {got} differ from the decided table {EXPECTED_STATUS_COUNTS}")
    for lib in NEPC_SUBTYPE_TERM:
        if lib not in fam.index:
            errs.append(f"{lib} (NEPC subtype library) is not in {MCRPC_STATES}")
    for lib, want in EXPECTED_DISEASE.items():
        if lib in per.index and per.loc[lib, "disease_ontology_term_id"] != want:
            errs.append(f"{lib} disease {per.loc[lib, 'disease_ontology_term_id']} (decided {want})")
    n_sub = int(((sm["disease_state"] == "primary_tumor") & (sm["disease_subtype"] != "")).sum())
    if n_sub != N_PRIMARY_SUBTYPE_LIBRARIES:
        errs.append(f"{n_sub} primary_tumor libraries with a disease_subtype (plan: {N_PRIMARY_SUBTYPE_LIBRARIES})")
    for lib in DEV_STAGE_BY_LIBRARY:
        if lib not in per.index:
            errs.append(f"{lib} (decided development stage) not in sample_metadata.csv")
    for lib in list(TISSUE_BY_LIBRARY) + list(NEPC_SUBTYPE_TERM) + list(ASSAY_BY_LIBRARY):
        if lib not in per.index:
            errs.append(f"{lib} (decided term) not in sample_metadata.csv")
    # by-library lists against the library's own sample_metadata.csv values
    for lib, (t, _) in NEPC_SUBTYPE_TERM.items():
        if lib in sm.index and sm.loc[lib, "neuroendocrine_type"] != NEPC_EXPECTED_NEUROENDOCRINE_TYPE[t]:
            errs.append(f"{lib} neuroendocrine_type '{sm.loc[lib, 'neuroendocrine_type']}' "
                        f"(the {t} list expects '{NEPC_EXPECTED_NEUROENDOCRINE_TYPE[t]}')")
    ne_values = set(NEPC_EXPECTED_NEUROENDOCRINE_TYPE.values())
    extra = sorted(set(sm.index[sm["neuroendocrine_type"].isin(ne_values)]) - set(NEPC_SUBTYPE_TERM))
    if extra:
        errs.append(f"libraries with neuroendocrine_type in {sorted(ne_values)} but no NEPC subtype term: {extra}")
    wong = set(sm.index[sm["study_id"] == WONG_STUDY_ID])
    if wong != set(DEV_STAGE_BY_LIBRARY):
        errs.append(f"libraries of study {WONG_STUDY_ID} {sorted(wong)} differ from the decided development-stage "
                    f"list {sorted(DEV_STAGE_BY_LIBRARY)}")
    if errs:
        stop("decided per-library terms not met:\n  " + "\n  ".join(errs))


# ----------------------------------------------------------------------------
# Per-library comparison with the earlier sample_level_schema_fields.csv
# (development stage added). Any other difference stops S5 and is listed.
# ----------------------------------------------------------------------------
def _libs(*parts):
    out = []
    for p in parts:
        if isinstance(p, tuple):
            out += [f"L{i:03d}" for i in range(p[0], p[1] + 1)]
        else:
            out.append(p)
    return set(out)


V2_EXPECTED_DIFFERENCES = {
    "assay_ontology_term_id": _libs((1, 10), "L083", (119, 125), "L132", "L133", (135, 139)),            # 25
    "disease_ontology_term_id": _libs("L078", "L079", "L082", "L083", (101, 109), "L110", "L135", "L136"),  # 16
    "development_stage_ontology_term_id": _libs("L128", "L129"),                                         # 2
    # L146 and L151 (UBERON:8410026 by library) are not on this list: if the earlier
    # sample_level_schema_fields.csv gives them UBERON:0002367, compare_with_v2 stops step 2 (a safe stop).
    # Check the earlier values against the decided terms and change this list only from that source.
    "tissue_ontology_term_id": _libs("L067", "L068", (101, 109)),                                        # 11
    "self_reported_ethnicity_ontology_term_id": _libs("L138"),                                           # 1
}
V2_EXPECTED_DONOR_ID_DIFFERENCES = 8  # the Hirz HP0-4 and P21 renames; the plan does not name the libraries
V2_EXPECTED_NEW = {"L071", "L126"}
V2_EXPECTED_GONE = {"L134", "L161"}

# ----------------------------------------------------------------------------
# Public-text scan
# applied to the text this code writes into the mapping tables and uns['cellxgene_conversion'].
# ----------------------------------------------------------------------------
BRITISH_PATTERNS = [
    r"tumour", r"labell(ed|ing)", r"normalis", r"colour", r"artefact", r"harmonis", r"\banalys(e|ed|ing)\b",
    r"\bcentre", r"licence", r"modell(ed|ing)", r"behaviour", r"programme", r"catalogue",
    r"organis(e|ed|ation)", r"characteris(e|ed|ation)", r"recognis(e|ed)", r"summaris(e|ed)",
    r"\bgrey\b", r"haemat", r"oestrogen", r"foetal", r"paediatric", r"ageing",
]
INTERNAL_PATTERNS = [
    r"/mctp", r"/scratch", r"/tmp", r"/home/", r"repository_files_v", r"cellxgene_schema_v", r"release_build",
    r"\.py\b", r"\bRB-\d+", r"\bH-\d+\b", r"(?<![A-Za-z])SVMs?(?![A-Za-z])",
    # the internal mapping and check tables of the run (none is deposited)
    r"\b(?:cell_type_ontology_mapping|disease_ontology_mapping|tissue_ontology_mapping|assay_ontology_mapping|development_stage_mapping|self_reported_ethnicity_mapping|sample_level_schema_fields|per_library_term_check)\.csv\b",
    # a role name and internal planning or record documents
    r"\bcoordinator\b", r"PLAN_v\d", r"README_P\d", r"Q4_terms_to_confirm", r"release_obs_columns",
]


def public_text_problems(text):
    t = str(text)
    hits = [p for p in BRITISH_PATTERNS if re.search(p, t, flags=re.IGNORECASE)]
    hits += [p for p in INTERNAL_PATTERNS if re.search(p, t)]
    return hits


def scan_table(df, name, strict_columns):
    """Stop on a hit in a column this code writes (strict); list hits in verbatim source columns."""
    strict, report = [], []
    for c in df.columns:
        # every non-numeric column, whatever its dtype: object, string or categorical
        if pd.api.types.is_bool_dtype(df[c].dtype) or pd.api.types.is_numeric_dtype(df[c].dtype):
            continue
        for i, v in df[c].items():
            if v is None or (not isinstance(v, str) and pd.isna(v)):
                continue
            h = public_text_problems(v)
            if h:
                (strict if c in strict_columns else report).append(f"{name} row {i} column {c}: {h}")
    for line in report:
        print(f"public-text note (verbatim source text, listed for review): {line}")
    if strict:
        stop("public-text scan failed:\n  " + "\n  ".join(strict))


# ----------------------------------------------------------------------------
# Deposit tables and the final hand-over. The objects are
# prepared for the Zenodo deposit only; nothing is submitted to the CELLxGENE portal.
# ----------------------------------------------------------------------------
# No table that build_mapping_tables.py (not included) writes is deposited. All stay internal records of the run.
DEPOSIT_MAPPING_TABLES = ()
AUTHOR_REVIEW_COLUMN = "author_review_needed"
REVIEW_FLAG_TABLES = ("cell_type_ontology_mapping.csv", "disease_ontology_mapping.csv", "tissue_ontology_mapping.csv")
# The author confirms the flagged terms,
# and the author_review_needed column is then left out of every table. PENDING until
# the confirmation, then a record of it. Until then the tables keep the column, and
# check_handover() blocks the final hand-over while the column says 'yes'.
AUTHOR_REVIEW_CONFIRMED = "confirmed by the author"
# "<record>; <md5 of the reviewed count>" of the reviewed count, in the S4 objects, of the first-version
# 'Neutrophils' cells that keep the S4 'Monocyte' label. PENDING until the count is reviewed. Monocyte's
# flag stays 'yes' until then. Afterwards it is set to 'no' from that count, or the column is dropped after the author's confirmation.
MONOCYTE_S4_COUNT_FROM = ("reviewed S4 count; caeb8f8120997d9b1307539458feb5fe; first-submission scope: of the 498 cells with the "
                          "first-submission label 'Neutrophils', 494 are released, 0 with the release cell_type "
                          "'Monocyte' (all 494 'Neutrophil'), and 4 are excluded")


def check_review_flags():
    """Consistency of the hand-over values: a set hand-over value has its form; Monocyte's flag stays 'yes' and the column
    stays in the tables while the reviewed S4 Monocyte count is PENDING. Runs in preflight(), build_mapping_tables.py (not included)
    and check_handover(); a PENDING value alone does not stop here."""
    bad = []
    for n in HANDOVER_NEEDED:
        v = globals()[n]
        if not _is_pending(v) and not DECIDED_FORMS[n][0](v):
            bad.append(f"{n} must be {DECIDED_FORMS[n][1]}")
    if _is_pending(MONOCYTE_S4_COUNT_FROM):
        if CELL_TYPE_MAP["Monocyte"][3] != "yes":
            bad.append("Monocyte's author_review_needed stays 'yes' while MONOCYTE_S4_COUNT_FROM is PENDING")
        if not _is_pending(AUTHOR_REVIEW_CONFIRMED):
            bad.append("AUTHOR_REVIEW_CONFIRMED is set while MONOCYTE_S4_COUNT_FROM is PENDING: the column is dropped "
                       "only after the reviewed Monocyte count of the source object")
    if bad:
        stop("hand-over values: " + "; ".join(bad))


def drop_author_review_column(tables):
    """The tables (name -> DataFrame) without the author_review_needed column once the author has confirmed the
    flagged terms (AUTHOR_REVIEW_CONFIRMED set); unchanged while it is PENDING."""
    check_review_flags()
    if _is_pending(AUTHOR_REVIEW_CONFIRMED):
        return tables
    return {n: (df.drop(columns=[AUTHOR_REVIEW_COLUMN]) if AUTHOR_REVIEW_COLUMN in df.columns else df)
            for n, df in tables.items()}


def check_handover(out_dir):
    """Final hand-over guard: reads the tables of a run folder and writes nothing. Stops while a
    hand-over value is PENDING, a deposit table is missing (none is listed), a review-flag table
    (REVIEW_FLAG_TABLES) is missing, any table's author_review_needed column says 'yes', or a
    table still has that column after the author's confirmation (built before it: rebuild the tables)."""
    require_decided(PIPELINE_NEEDED)
    check_forms(PIPELINE_NEEDED)
    check_review_flags()
    problems = []
    pend = [n for n in HANDOVER_NEEDED if _is_pending(globals()[n])]
    if pend:
        problems.append(f"PENDING {pend} (the author's confirmation of the flagged terms; the reviewed Monocyte "
                        "count of the source object)")
    for name in DEPOSIT_MAPPING_TABLES:
        if not os.path.isfile(os.path.join(out_dir, name)):
            problems.append(f"deposit table {name} is not in {out_dir}")
    for name in sorted(set(REVIEW_FLAG_TABLES) | set(DEPOSIT_MAPPING_TABLES)):
        path = os.path.join(out_dir, name)
        if not os.path.isfile(path):
            # its flags cannot be read, so a missing review-flag table stops
            if name in REVIEW_FLAG_TABLES:
                problems.append(f"review-flag table {name} is not in {out_dir}")
            continue
        df = pd.read_csv(path, dtype=str, keep_default_na=False)
        if AUTHOR_REVIEW_COLUMN not in df.columns:
            continue
        yes = df[df[AUTHOR_REVIEW_COLUMN] == "yes"]
        if len(yes):
            problems.append(f"{name}: {AUTHOR_REVIEW_COLUMN} says 'yes' for {sorted(yes.iloc[:, 0].astype(str))}")
        elif not _is_pending(AUTHOR_REVIEW_CONFIRMED):
            problems.append(f"{name} still has the {AUTHOR_REVIEW_COLUMN} column after the confirmation "
                            f"({AUTHOR_REVIEW_CONFIRMED}); rebuild the tables")
    if problems:
        stop("final hand-over blocked:\n  " + "\n  ".join(problems))
    print(f"hand-over check passed for {out_dir}: deposit tables {list(DEPOSIT_MAPPING_TABLES)} and review-flag tables "
          f"{list(REVIEW_FLAG_TABLES)} present; no "
          f"{AUTHOR_REVIEW_COLUMN} column; {HANDOVER_NEEDED} recorded", flush=True)
