# Code release header: tests pairs of GEM wells within each study for excess sharing of Cell Ranger-called cell barcodes. Three studies of one lab form one group.
# Input: a sample manifest (placeholder) and the filtered barcode list of each GEM well's Cell Ranger output.
# Output: barcode_sharing_pairs.csv and barcode_sharing_summary.csv, in the script's folder.
# Expected shared barcodes = nA x nB / W, with W the 10x whitelist size of the chemistry. The p value is the hypergeometric upper tail.
# Excess sharing: at least 5-fold over expected and p < 1e-6. Pairs with different whitelists are not tested. Methods: Raw sequencing processing.
# In the run, pairs met both cutoffs only in the four barcode-swap studies. No tested pair of the other studies exceeded 4.7-fold.
"""Excess sharing of Cell Ranger-called cell barcodes between libraries (GEM wells) sequenced in the same study
(the Strand lab studies P15/P16/P17 are tested together). Expected shared barcodes under independent draws = nA*nB/W (W = size of the 10x
barcode whitelist of the chemistry Cell Ranger detected: 737,280 for 737K-august-2016 (3' v2, 5' v1/PE), 6,794,880 for 3M-february-2018 (3' v3);
the v3 size is the 10x-documented value, not counted from a local file), hypergeometric upper-tail p, as for MM033 (3,298 shared vs ~16 expected).
Pairs with different whitelists are listed but not tested. Read-only on all Cell Ranger outputs; writes tables here."""
import csv, glob, gzip, itertools, os, collections
import h5py
from scipy.stats import hypergeom
R = "<PROJECT_DIR>"
OUT = os.path.dirname(os.path.abspath(__file__))
rows = list(csv.DictReader(open("<input: library manifest, one row per library (columns lib, study, chemistry_by_data, redo_input_dirs)>")))
override = {"L111": "<CELLRANGER_RERUN_DIR>/Crowley/MM033_allflowcells",
            "L126": "<CELLRANGER_RERUN_DIR>/Henry_GSE117403/D17PrPzF_Via_rerun"}
units = collections.OrderedDict()
for r in rows:
    dirs = [override[r["lib"]]] if r["lib"] in override else [x.strip() for x in r["redo_input_dirs"].split(";")]
    grp = "P15_P16_P17" if r["study"][:3] in ("P15", "P16", "P17") else r["study"]
    for d in dirs:  # Tabula Sapiens: one unit per 10x channel
        u = units.setdefault(d, {"libs": [], "group": grp})
        u["chem"] = r["chemistry_by_data"].split(" ")[0]
        u["libs"].append(r["lib"] + ("_" + d.rsplit("_", 2)[-2] + "_" + d.rsplit("_", 1)[-1] if len(dirs) > 1 else ""))
for d, u in units.items():
    mi = (glob.glob(d + "/outs/molecule_info.h5") + glob.glob(d + "/outs/per_sample_outs/*/count/sample_molecule_info.h5"))[0]
    fb = (glob.glob(d + "/outs/filtered_feature_bc_matrix/barcodes.tsv.gz") + glob.glob(d + "/outs/per_sample_outs/*/count/sample_filtered_feature_bc_matrix/barcodes.tsv.gz"))[0]
    u["W"] = {"SC3Pv2": 737280, "SC5P-PE": 737280, "SC3Pv3": 6794880}[u["chem"]]
    u["bc"] = set(l.strip()[:16] for l in gzip.open(fb, "rt"))
    u["name"] = "+".join(u["libs"])
out = []
by = collections.defaultdict(list)
for d, u in units.items():
    by[u["group"]].append(u)
for g, us in by.items():
    for a, b in itertools.combinations(us, 2):
        s = len(a["bc"] & b["bc"])
        if a["W"] != b["W"]:
            out.append(dict(group=g, lib_a=a["name"], lib_b=b["name"], n_a=len(a["bc"]), n_b=len(b["bc"]), whitelist_a=a["W"], whitelist_b=b["W"], shared=s, expected="", fold="", p_hypergeom="", note="different whitelists, not tested"))
            continue
        e = len(a["bc"]) * len(b["bc"]) / a["W"]
        p = hypergeom.sf(s - 1, a["W"], len(a["bc"]), len(b["bc"])) if s > 0 else 1.0
        out.append(dict(group=g, lib_a=a["name"], lib_b=b["name"], n_a=len(a["bc"]), n_b=len(b["bc"]), whitelist_a=a["W"], whitelist_b=b["W"], shared=s, expected=round(e, 2), fold=round(s / e, 1) if e else "", p_hypergeom="%.3g" % p, note=""))
F = list(out[0].keys())
tmp = OUT + "/barcode_sharing_pairs.csv.tmp"
with open(tmp, "w", newline="") as h:
    w = csv.DictWriter(h, fieldnames=F); w.writeheader(); w.writerows(out)
os.replace(tmp, OUT + "/barcode_sharing_pairs.csv")
summ = []
for g, us in by.items():
    pr = [x for x in out if x["group"] == g and x["note"] == ""]
    hits = [x for x in pr if x["fold"] != "" and float(x["fold"]) >= 5 and float(x["p_hypergeom"]) < 1e-6]
    summ.append(dict(group=g, units=len(us), pairs_tested=len(pr), pairs_untested=len([x for x in out if x["group"] == g and x["note"]]),
                     pairs_fold5_p1e6=len(hits), max_fold=max([float(x["fold"]) for x in pr if x["fold"] != ""], default=""),
                     max_shared=max([x["shared"] for x in pr], default=0)))
with open(OUT + "/barcode_sharing_summary.csv", "w", newline="") as h:
    w = csv.DictWriter(h, fieldnames=list(summ[0].keys())); w.writeheader(); w.writerows(summ)
for s in summ: print(s)
