"""Fase 4 — Spatial enrichment: alpha genetic signal in the target ROIs.

Tests whether the cortical alpha-source ROIs (Tabarelli 2022) carry a stronger
"alpha transcriptomic signal" than the rest of cortex, using a spin-test null
that preserves spatial autocorrelation.

Two per-region alpha scores (see enrichment_lib):
  primary  'continuous' — Spearman corr between a region's expression profile and
           the genome-wide gene-level Z (all shared genes; no p-value threshold).
  sensit.  'top100'     — mean z-scored expression of the top-N MAGMA gene set.
Two ROI sets: 'full' (41 Tabarelli ROIs = 82 parcels) and 'hub' (8C/8Av/9a).

Shared machinery (centroids, sphere spins, scores, spin p-value) lives in
enrichment_lib.py. Fase 5 (05_sensitivity.py) reuses it across all phenotypes.

Run:  python scripts/04_spin_test_enrichment.py --pheno alphaCz
Outputs: results/tables/enrichment_<pheno>.csv, results/figures/enrichment_<pheno>.png
"""
from __future__ import annotations

import argparse

import numpy as np
import pandas as pd

from utils import ROOT, load_config, log_run
import enrichment_lib as EL

RESULTS = ROOT / "results"


def main() -> None:
    cfg = load_config()
    ap = argparse.ArgumentParser()
    ap.add_argument("--pheno", default=cfg["gwas"]["primary"])
    ap.add_argument("--n-spins", type=int, default=cfg["enrichment"]["n_spins"])
    args = ap.parse_args()
    seed = cfg["seed"]
    alt = cfg["enrichment"]["alternative"]

    spins, ids = EL.load_or_make_spins(args.n_spins, seed)
    masks = EL.masks_for(ids)
    scores = EL.score_maps(args.pheno, ids, cfg["enrichment"]["corr"])

    rows = []
    for sname, svec in scores.items():
        for rname, mask in masks.items():
            res = EL.spin_pvalue(svec, mask, spins, alt)
            res.update({"phenotype": args.pheno, "score": sname, "roi_set": rname,
                        "n_target": int(mask.sum())})
            rows.append(res)
            print(f"  {sname:10s} {rname:4s}: diff={res['diff']:+.4f} "
                  f"z={res['z_vs_null']:+.2f} p_spin={res['p_spin']:.4f}")

    out = pd.DataFrame(rows)[["phenotype", "score", "roi_set", "n_target",
                              "mean_target", "mean_rest", "diff", "null_mean",
                              "null_sd", "z_vs_null", "p_spin"]]
    out_path = RESULTS / "tables" / f"enrichment_{args.pheno}.csv"
    out.to_csv(out_path, index=False)
    print(f"[write] {out_path}")

    _figure(scores, masks, spins, args.pheno, alt)
    log_run("04_spin_test_enrichment", {
        "phenotype": args.pheno, "n_spins": int(args.n_spins), "seed": seed,
        "alternative": alt, "corr": cfg["enrichment"]["corr"],
        "results": out.to_dict("records")})


def _figure(scores, masks, spins, pheno, alt) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    combos = [(s, r) for s in scores for r in masks]
    fig, axes = plt.subplots(len(scores), len(masks), figsize=(9, 7), squeeze=False)
    for ax, (sname, rname) in zip(axes.ravel(), combos):
        svec, mask = scores[sname], masks[rname]
        obs = svec[mask].mean()
        null = np.array([svec[spins[:, s]][mask].mean() for s in range(spins.shape[1])])
        ax.hist(null, bins=60, color="0.7", edgecolor="none")
        ax.axvline(obs, color="crimson", lw=2, label=f"observed={obs:.3f}")
        p = (1 + np.sum(null >= obs)) / (len(null) + 1) if alt == "greater" else np.nan
        ax.set_title(f"{sname} | {rname} ROIs  (p_spin={p:.4f})", fontsize=10)
        ax.set_xlabel("mean alpha score in target (spin null)")
        ax.legend(fontsize=8)
    fig.suptitle(f"Spatial enrichment of alpha genetic signal — {pheno}", fontsize=12)
    fig.tight_layout()
    path = RESULTS / "figures" / f"enrichment_{pheno}.png"
    fig.savefig(path, dpi=130)
    print(f"[write] {path}")


if __name__ == "__main__":
    main()
