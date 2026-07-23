"""Fase 5 — Specificity & sensitivity of the alpha enrichment.

Runs the Fase 4 spin-test enrichment across ALL phenotypes and assembles a single
comparison table, so the alpha result can be judged against:
  - a replication phenotype (alphaOcc), and
  - specificity controls (thetaCz, betaCz, deltaCz) that should be null.

Applies FDR correction across the full phenotype x score x ROI-set grid, and
flags whether alpha (and only alpha) is enriched.

Requires each phenotype's Fase 1 output (<pheno>_geneZ.csv, <pheno>_genes_top100.txt).

Run:  python scripts/05_sensitivity.py
Outputs: results/tables/enrichment_all_phenotypes.csv
         results/figures/specificity_summary.png
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

from utils import ROOT, load_config, log_run
import enrichment_lib as EL

DERIVED = ROOT / "data" / "derived"
RESULTS = ROOT / "results"


def available_phenos(cfg) -> list[str]:
    order = [cfg["gwas"]["primary"], "alphaOcc", "peakOcc", *cfg["gwas"]["controls"]]
    seen, out = set(), []
    for p in order:
        if p in seen:
            continue
        seen.add(p)
        if (DERIVED / f"{p}_geneZ.csv").exists():
            out.append(p)
        else:
            print(f"[skip] {p}: no Fase 1 output yet")
    return out


def main() -> None:
    cfg = load_config()
    n_spins = cfg["enrichment"]["n_spins"]
    seed = cfg["seed"]
    corr = cfg["enrichment"]["corr"]
    alt = cfg["enrichment"]["alternative"]

    phenos = available_phenos(cfg)
    if not phenos:
        raise SystemExit("No phenotype Fase 1 outputs found. Run scripts/01_* first.")

    spins, ids = EL.load_or_make_spins(n_spins, seed)
    masks = EL.masks_for(ids)

    band_of = {"alphaCz": "alpha-power", "alphaOcc": "alpha-power", "peakOcc": "alpha-frequency"}
    rows = []
    for ph in phenos:
        smaps = EL.score_maps(ph, ids, corr)
        for sname, svec in smaps.items():
            for rname, mask in masks.items():
                res = EL.spin_pvalue(svec, mask, spins, alt)
                res.update({"phenotype": ph, "band": band_of.get(ph, "control"),
                            "score": sname, "roi_set": rname, "n_target": int(mask.sum())})
                rows.append(res)

    df = pd.DataFrame(rows)
    df["p_fdr"] = multipletests(df["p_spin"], method="fdr_bh")[1]
    df = df[["phenotype", "band", "score", "roi_set", "n_target", "diff",
             "z_vs_null", "p_spin", "p_fdr"]].sort_values(
        ["score", "roi_set", "phenotype"]).reset_index(drop=True)

    out = RESULTS / "tables" / "enrichment_all_phenotypes.csv"
    df.to_csv(out, index=False)
    print(df.to_string(index=False))
    print(f"\n[write] {out}")

    _figure(df)
    log_run("05_sensitivity", {"phenotypes": phenos, "n_spins": int(n_spins),
                               "seed": seed, "results": df.to_dict("records")})


def _figure(df: pd.DataFrame) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sub = df[df.roi_set == "full"].copy()
    scores = sub.score.unique()
    fig, axes = plt.subplots(1, len(scores), figsize=(11, 4.5), squeeze=False)
    for ax, sc in zip(axes.ravel(), scores):
        d = sub[sub.score == sc].sort_values("phenotype")
        colors = ["#c0392b" if b == "alpha" else "#7f8c8d" for b in d.band]
        ax.bar(d.phenotype, -np.log10(d.p_spin), color=colors)
        ax.axhline(-np.log10(0.05), ls="--", c="k", lw=1, label="p=0.05")
        ax.set_title(f"{sc} score, 41-ROI set")
        ax.set_ylabel("-log10 p_spin")
        ax.tick_params(axis="x", rotation=45)
        ax.legend(fontsize=8)
    fig.suptitle("Specificity: alpha (red) vs control bands (grey)", fontsize=12)
    fig.tight_layout()
    path = RESULTS / "figures" / "specificity_summary.png"
    fig.savefig(path, dpi=130)
    print(f"[write] {path}")


if __name__ == "__main__":
    main()
