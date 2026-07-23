"""Method-validation analyses for the spatial enrichment test.

(A) Positive control: show the spin/enrichment machinery DOES detect a true
    spatial effect. We use the dominant AHBA expression gradient (PC1) as a
    known-strong signal: with PC1 gene loadings as the gene statistic and the
    top-tertile PC1 regions as the target, the continuous score must recover a
    highly significant enrichment. If it does, the alpha null is not a floor
    effect of the method.

(B) Gene-set (co-expression-aware) null: beyond the spin (which permutes space),
    test whether the alpha top-100 set is special vs random gene sets of equal
    size. p_geneset = fraction of random sets whose target-vs-rest difference
    >= observed. Addresses Fulcher et al. (2021).

(C) Spin validation: compare the custom sphere spin against brainsmash
    spatial-autocorrelation-preserving surrogate maps (Euclidean parcel
    distances on the fsaverage5 white surface) for the primary alphaCz
    continuous map. Concordant p-values validate the custom spin.

Run:  python scripts/07_method_validation.py
Outputs: results/tables/validation_positive_control.csv,
         results/tables/validation_geneset_null.csv,
         results/tables/validation_spin_validation.csv
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import nibabel as nib

from utils import ROOT, load_config, log_run
import enrichment_lib as EL

DERIVED = ROOT / "data" / "derived"
REFERENCE = ROOT / "data" / "reference"
RESULTS = ROOT / "results"


def load_common():
    cfg = load_config()
    spins, ids = EL.load_or_make_spins(cfg["enrichment"]["n_spins"], cfg["seed"])
    masks = EL.masks_for(ids)
    expr = pd.read_csv(DERIVED / "ahba_expression_glasser.csv", index_col=0).loc[ids]
    return cfg, spins, ids, masks, expr


# ---------------------------- (A) positive control ----------------------------

def positive_control(cfg, spins, ids, expr) -> pd.DataFrame:
    E = expr.to_numpy()
    Ec = E - E.mean(0, keepdims=True)
    U, S, Vt = np.linalg.svd(Ec, full_matrices=False)
    pc1_region = U[:, 0] * S[0]                    # region scores on PC1
    pc1_loadings = Vt[0, :]                         # gene loadings
    # orient so the target end is positive
    if np.corrcoef(pc1_region, (E @ pc1_loadings))[0, 1] < 0:
        pc1_region, pc1_loadings = -pc1_region, -pc1_loadings
    # target = top-tertile PC1 regions (a genuine, strong spatial signal)
    thr = np.quantile(pc1_region, 2 / 3)
    target = pc1_region >= thr
    geneZ = pd.DataFrame({"symbol": expr.columns, "ZSTAT": pc1_loadings, "P": np.nan})
    score = EL.continuous_score(expr, geneZ, cfg["enrichment"]["corr"])
    res = EL.spin_pvalue(np.asarray(score), target, spins, "greater")
    res.update({"analysis": "positive_control_PC1", "n_target": int(target.sum())})
    return pd.DataFrame([res])[["analysis", "n_target", "diff", "z_vs_null", "p_spin"]]


# ---------------------------- (B) gene-set null ----------------------------

def geneset_null(cfg, ids, masks, expr, n_null=10000) -> pd.DataFrame:
    rng = np.random.default_rng(cfg["seed"])
    target = masks["full"]
    genes = expr.columns.to_numpy()
    alpha = [g for g in (DERIVED / "alphaCz_genes_top100.txt").read_text().split() if g in expr.columns]
    def diff_for(gene_list):
        s = EL.topn_score(expr, list(gene_list))
        return s[target].mean() - s[~target].mean()
    obs = diff_for(alpha)
    k = len(alpha)
    null = np.empty(n_null)
    idx = np.arange(len(genes))
    for i in range(n_null):
        null[i] = diff_for(genes[rng.choice(idx, size=k, replace=False)])
    p = (1 + np.sum(null >= obs)) / (n_null + 1)
    res = {"analysis": "geneset_null_top100_full", "n_genes": k, "obs_diff": float(obs),
           "null_mean": float(null.mean()), "null_sd": float(null.std()),
           "z_vs_null": float((obs - null.mean()) / null.std()), "p_geneset": float(p)}
    return pd.DataFrame([res])


# ---------------------------- (C) spin validation ----------------------------

def white_centroids(ids):
    import abagen
    fs = abagen.images.fetch_fsaverage5()
    labs = {"lh": np.asarray(nib.load(REFERENCE / "glasser_fsaverage5_lh.label.gii").agg_data()).astype(int),
            "rh": np.asarray(nib.load(REFERENCE / "glasser_fsaverage5_rh.label.gii").agg_data()).astype(int)}
    verts = {"lh": np.asarray(fs.lh.vertices), "rh": np.asarray(fs.rh.vertices)}
    cent = {}
    for h in ("lh", "rh"):
        for i in np.unique(labs[h][labs[h] != 0]):
            cent[int(i)] = verts[h][labs[h] == i].mean(0)
    return np.array([cent[i] for i in ids])


def spin_validation(cfg, spins, ids, masks, expr, n_surr=1000) -> pd.DataFrame:
    from scipy.spatial.distance import cdist
    from brainsmash.mapgen.base import Base
    geneZ = pd.read_csv(DERIVED / "alphaCz_geneZ.csv")
    score = np.asarray(EL.continuous_score(expr, geneZ, cfg["enrichment"]["corr"]))
    target = masks["full"]
    obs = score[target].mean() - score[~target].mean()
    # custom spin p (reference)
    spin_res = EL.spin_pvalue(score, target, spins, "greater")
    # brainsmash surrogates on Euclidean white-surface parcel distances
    D = cdist(white_centroids(ids), white_centroids(ids))
    surr = Base(x=score, D=D, seed=cfg["seed"])(n=n_surr)
    null = np.array([s[target].mean() - s[~target].mean() for s in surr])
    p_bs = (1 + np.sum(null >= obs)) / (n_surr + 1)
    res = {"analysis": "spin_validation_alphaCz_continuous", "obs_diff": float(obs),
           "p_spin_custom": float(spin_res["p_spin"]),
           "p_brainsmash": float(p_bs), "n_brainsmash_surrogates": n_surr}
    return pd.DataFrame([res])


def main() -> None:
    cfg, spins, ids, masks, expr = load_common()

    pc = positive_control(cfg, spins, ids, expr)
    pc.to_csv(RESULTS / "tables" / "validation_positive_control.csv", index=False)
    print("[A] positive control\n", pc.to_string(index=False))

    gs = geneset_null(cfg, ids, masks, expr)
    gs.to_csv(RESULTS / "tables" / "validation_geneset_null.csv", index=False)
    print("\n[B] gene-set null\n", gs.to_string(index=False))

    sv = spin_validation(cfg, spins, ids, masks, expr)
    sv.to_csv(RESULTS / "tables" / "validation_spin_validation.csv", index=False)
    print("\n[C] spin validation\n", sv.to_string(index=False))

    log_run("07_method_validation", {
        "positive_control": pc.to_dict("records"),
        "geneset_null": gs.to_dict("records"),
        "spin_validation": sv.to_dict("records")})


if __name__ == "__main__":
    main()
