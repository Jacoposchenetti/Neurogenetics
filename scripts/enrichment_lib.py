"""Shared enrichment machinery for Fase 4 (04_spin_test_enrichment) and Fase 5
(05_sensitivity): parcel centroids, sphere spin nulls, alpha scores, spin p-value.

Kept as a plain-named module so both entry points can import it (script files
starting with a digit are not importable).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import nibabel as nib
from scipy.spatial import cKDTree
from scipy.stats import rankdata

from netneurotools.datasets import fetch_hcp_standards

from utils import ROOT

DERIVED = ROOT / "data" / "derived"
REFERENCE = ROOT / "data" / "reference"
RESULTS = ROOT / "results"

# 41 Tabarelli ROIs (labels), mirrored from scripts/02_roi_mapping.py
ROI_LABELS = [
    "V1", "V2", "ProS", "V3", "V4", "V6", "V6A", "V7", "IPS1", "V3A", "V3B",
    "V3CD", "IP0", "PGp", "LO1", "LO2", "1", "2", "3a", "3b", "4", "6mp", "6d",
    "8BL", "9p", "9m", "9a", "8Ad", "9-46d", "8BM", "8Av", "46", "8C", "p9-46v",
    "a32pr", "d32", "a9-46v", "10d", "p10p", "p47r", "IFSa",
]
HUB_LABELS = ["8C", "8Av", "9a"]


def parcel_centroids() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(ids, centroids_xyz, hemi) for the 360 Glasser parcels on the fsaverage5 sphere.

    ids sorted ascending (1..360); hemi is 'L'/'R' per parcel.
    """
    std = Path(fetch_hcp_standards()) / "resample_fsaverage"
    spheres = {"lh": std / "fsaverage5_std_sphere.L.10k_fsavg_L.surf.gii",
               "rh": std / "fsaverage5_std_sphere.R.10k_fsavg_R.surf.gii"}
    labels_f = {"lh": REFERENCE / "glasser_fsaverage5_lh.label.gii",
                "rh": REFERENCE / "glasser_fsaverage5_rh.label.gii"}
    ids, cents, hemi = [], [], []
    for h, hlab in (("lh", "L"), ("rh", "R")):
        xyz = np.asarray(nib.load(str(spheres[h])).agg_data("NIFTI_INTENT_POINTSET"))
        lab = np.asarray(nib.load(str(labels_f[h])).agg_data()).astype(int)
        for i in np.unique(lab[lab != 0]):
            c = xyz[lab == i].mean(0)
            ids.append(int(i)); cents.append(c / np.linalg.norm(c)); hemi.append(hlab)
    order = np.argsort(ids)
    return np.asarray(ids)[order], np.asarray(cents)[order], np.asarray(hemi)[order]


def _rand_rotation(rng) -> np.ndarray:
    q, r = np.linalg.qr(rng.normal(size=(3, 3)))
    q = q @ np.diag(np.sign(np.diag(r)))
    if np.linalg.det(q) < 0:
        q[:, 0] = -q[:, 0]
    return q


def gen_spins(cents: np.ndarray, hemi: np.ndarray, n_spins: int, seed: int) -> np.ndarray:
    """(n_parcels, n_spins) permutation indices via sphere rotation (nearest-neighbour).

    Left hemisphere is rotated; right hemisphere gets the x-reflected rotation so
    the two hemispheres stay mirror-consistent (Alexander-Bloch 2018).
    """
    rng = np.random.default_rng(seed)
    L = np.where(hemi == "L")[0]
    R = np.where(hemi == "R")[0]
    refl = np.diag([-1.0, 1.0, 1.0])
    treeL, treeR = cKDTree(cents[L]), cKDTree(cents[R])
    spins = np.empty((len(cents), n_spins), dtype=np.int32)
    for s in range(n_spins):
        rot = _rand_rotation(rng)
        spins[L, s] = L[treeL.query(cents[L] @ rot)[1]]
        spins[R, s] = R[treeR.query((cents[R] @ refl) @ rot @ refl)[1]]
    return spins


def load_or_make_spins(n_spins: int, seed: int) -> tuple[np.ndarray, np.ndarray]:
    """Return (spins, ids) using a cached .npy when available."""
    ids, cents, hemi = parcel_centroids()
    path = RESULTS / "derived" / f"spins_glasser_fsaverage5_{n_spins}_s{seed}.npy"
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        spins = np.load(path)
    else:
        spins = gen_spins(cents, hemi, n_spins, seed)
        np.save(path, spins)
    return spins, ids


def continuous_score(expr: pd.DataFrame, geneZ: pd.DataFrame, method: str) -> np.ndarray:
    """Per-region association between expression profile and gene-level Z."""
    z = geneZ.dropna(subset=["ZSTAT"]).drop_duplicates("symbol").set_index("symbol")["ZSTAT"]
    shared = [g for g in expr.columns if g in z.index]
    E = expr[shared].to_numpy()
    zz = z.loc[shared].to_numpy()
    if method == "spearman":
        E = np.apply_along_axis(rankdata, 1, E)
        zz = rankdata(zz)
    E = E - E.mean(1, keepdims=True)
    zz = zz - zz.mean()
    return (E @ zz) / np.sqrt((E ** 2).sum(1) * (zz ** 2).sum())


def topn_score(expr: pd.DataFrame, genes: list[str]) -> np.ndarray:
    present = [g for g in genes if g in expr.columns]
    sub = expr[present]
    z = (sub - sub.mean(0)) / sub.std(0, ddof=0)
    return z.mean(1).to_numpy()


def spin_pvalue(score: np.ndarray, mask: np.ndarray, spins: np.ndarray,
                alternative: str = "greater") -> dict:
    obs = score[mask].mean()
    rest = score[~mask].mean()
    null = np.array([score[spins[:, s]][mask].mean() for s in range(spins.shape[1])])
    if alternative == "greater":
        p = (1 + np.sum(null >= obs)) / (len(null) + 1)
    elif alternative == "less":
        p = (1 + np.sum(null <= obs)) / (len(null) + 1)
    else:
        p = (1 + np.sum(np.abs(null - null.mean()) >= abs(obs - null.mean()))) / (len(null) + 1)
    return {"mean_target": float(obs), "mean_rest": float(rest),
            "diff": float(obs - rest), "null_mean": float(null.mean()),
            "null_sd": float(null.std()),
            "z_vs_null": float((obs - null.mean()) / null.std()), "p_spin": float(p)}


def masks_for(ids: np.ndarray) -> dict:
    """Boolean target masks (aligned to sorted ids) for the ROI sets."""
    info = pd.read_csv(REFERENCE / "glasser_360_info.csv").set_index("id")["label"]
    labels = info.loc[ids].to_numpy()
    return {"full": np.isin(labels, ROI_LABELS), "hub": np.isin(labels, HUB_LABELS)}


def score_maps(pheno: str, ids: np.ndarray, corr: str) -> dict:
    """Both alpha score maps for a phenotype, aligned to sorted ids."""
    expr = pd.read_csv(DERIVED / "ahba_expression_glasser.csv", index_col=0).loc[ids]
    geneZ = pd.read_csv(DERIVED / f"{pheno}_geneZ.csv")
    top_genes = (DERIVED / f"{pheno}_genes_top100.txt").read_text().split()
    return {"continuous": continuous_score(expr, geneZ, corr),
            "top100": topn_score(expr, top_genes)}
