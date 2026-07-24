"""Fase 6 — Parcellation sensitivity of the alpha enrichment.

Re-runs the spin-test enrichment on alternative fsaverage5 parcellations to test
whether granularity changes the (null) result, motivated by the resolution
mismatch between MNE-reconstructed alpha sources (~cm point-spread) and fine
transcriptomic parcels:
  schaefer100  coarse (100 parcels) — closer to MNE spatial resolution
  yan600       fine (600 parcels)   — Yan2023 homotopic (Kong2022 17-networks)

The 41 Glasser alpha-source ROIs are re-mapped onto each parcellation by
fsaverage5 surface overlap (a parcel is 'target' if >= OVERLAP_THR of its
cortical vertices fall in the Glasser alpha-source set).

Run:  python scripts/06_parcellation_sensitivity.py --parc schaefer100 --pheno alphaCz
      python scripts/06_parcellation_sensitivity.py --parc yan600 --pheno alphaCz
Outputs: data/reference/<parc>_fsaverage5_{lh,rh}.label.gii, <parc>_info.csv
         data/derived/ahba_expression_<parc>.csv (cached; abagen ~10-15 min)
         appends to results/tables/enrichment_parcellations.csv
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import nibabel as nib
from netneurotools.datasets import fetch_hcp_standards, fetch_schaefer2018

from utils import ROOT, load_config, log_run
import enrichment_lib as EL

REFERENCE = ROOT / "data" / "reference"
DERIVED = ROOT / "data" / "derived"
RESULTS = ROOT / "results"
OVERLAP_THR = 0.5


YAN600_SRC = (
    "https://raw.githubusercontent.com/ThomasYeoLab/CBIG/master/stable_projects/"
    "brain_parcellation/Yan2023_homotopic/parcellations/FreeSurfer/fsaverage5/label/kong17"
)


def annot_paths(parc: str) -> tuple[str, str]:
    if parc == "schaefer100":
        s = fetch_schaefer2018(version="fsaverage5")["100Parcels7Networks"]
        return str(s[0]), str(s[1])
    if parc == "yan600":
        lh = REFERENCE / "yan600_lh.annot"
        rh = REFERENCE / "yan600_rh.annot"
        missing = [p.name for p in (lh, rh) if not p.exists()]
        if missing:
            raise SystemExit(
                f"Yan-600 parcellation not found in data/reference/ (missing: {', '.join(missing)}).\n"
                "It is not redistributed here. Download the two fsaverage5 annot files from CBIG\n"
                "and rename them (the target names differ from the source names):\n"
                f"  curl -L {YAN600_SRC}/lh.600Parcels_Kong2022_17Networks.annot \\\n"
                "       -o data/reference/yan600_lh.annot\n"
                f"  curl -L {YAN600_SRC}/rh.600Parcels_Kong2022_17Networks.annot \\\n"
                "       -o data/reference/yan600_rh.annot\n"
                "See data/reference/README.md for details."
            )
        return str(lh), str(rh)
    raise SystemExit(f"unknown parcellation: {parc}")


def global_labels(parc: str):
    """Per-vertex global integer labels (lh,rh) + atlas_info (id,label,hemisphere).

    LH parcels get ids 1..nL, RH parcels nL+1..nL+nR; 0 = medial wall/background.
    """
    lh, rh = annot_paths(parc)
    al = nib.freesurfer.read_annot(lh)
    ar = nib.freesurfer.read_annot(rh)
    lab_l, names_l = al[0].astype(int), al[2]
    lab_r, names_r = ar[0].astype(int), ar[2]
    nL = int(lab_l.max())
    gl = lab_l.copy()
    gr = np.where(lab_r > 0, lab_r + nL, 0)
    rows = []
    for k in range(1, nL + 1):
        nm = names_l[k].decode() if isinstance(names_l[k], bytes) else names_l[k]
        rows.append({"id": k, "label": nm, "hemisphere": "L", "structure": "cortex"})
    for k in range(1, int(lab_r.max()) + 1):
        nm = names_r[k].decode() if isinstance(names_r[k], bytes) else names_r[k]
        rows.append({"id": k + nL, "label": nm, "hemisphere": "R", "structure": "cortex"})
    return {"lh": gl, "rh": gr}, pd.DataFrame(rows)


def write_label_gifti(labels: np.ndarray, out: Path) -> None:
    darr = nib.gifti.GiftiDataArray(labels.astype(np.int32),
                                    intent="NIFTI_INTENT_LABEL", datatype="NIFTI_TYPE_INT32")
    nib.gifti.GiftiImage(darrays=[darr]).to_filename(str(out))


def alpha_source_vertices():
    """Boolean per-vertex mask (lh,rh) of Glasser alpha-source ROI membership."""
    info = pd.read_csv(REFERENCE / "glasser_360_info.csv")
    alpha_ids = set(info.loc[info["label"].isin(EL.ROI_LABELS), "id"])
    out = {}
    for h in ("lh", "rh"):
        lab = np.asarray(nib.load(REFERENCE / f"glasser_fsaverage5_{h}.label.gii").agg_data()).astype(int)
        out[h] = np.isin(lab, list(alpha_ids))
    return out


def parc_centroids_and_mask(glab: dict, info: pd.DataFrame):
    """Sphere centroids per parcel (sorted ids) + alpha-overlap target mask."""
    std = Path(fetch_hcp_standards()) / "resample_fsaverage"
    sph = {"lh": std / "fsaverage5_std_sphere.L.10k_fsavg_L.surf.gii",
           "rh": std / "fsaverage5_std_sphere.R.10k_fsavg_R.surf.gii"}
    averts = alpha_source_vertices()
    ids, cents, hemi, frac = [], [], [], []
    for h, hlab in (("lh", "L"), ("rh", "R")):
        xyz = np.asarray(nib.load(str(sph[h])).agg_data("NIFTI_INTENT_POINTSET"))
        lab = glab[h]
        for i in np.unique(lab[lab != 0]):
            v = lab == i
            c = xyz[v].mean(0)
            ids.append(int(i)); cents.append(c / np.linalg.norm(c)); hemi.append(hlab)
            frac.append(averts[h][v].mean())
    order = np.argsort(ids)
    ids = np.asarray(ids)[order]
    cents = np.asarray(cents)[order]
    hemi = np.asarray(hemi)[order]
    frac = np.asarray(frac)[order]
    return ids, cents, hemi, (frac >= OVERLAP_THR)


def build_and_run_abagen(parc: str, glab: dict, info: pd.DataFrame, cfg) -> Path:
    out = DERIVED / f"ahba_expression_{parc}.csv"
    if out.exists():
        print(f"[cache] {out}")
        return out
    lh_g = REFERENCE / f"{parc}_fsaverage5_lh.label.gii"
    rh_g = REFERENCE / f"{parc}_fsaverage5_rh.label.gii"
    write_label_gifti(glab["lh"], lh_g)
    write_label_gifti(glab["rh"], rh_g)
    info.to_csv(REFERENCE / f"{parc}_info.csv", index=False)
    import abagen
    a = cfg["abagen"]
    expr = abagen.get_expression_data(
        (str(lh_g), str(rh_g)), atlas_info=REFERENCE / f"{parc}_info.csv",
        ibf_threshold=a["ibf_threshold"], probe_selection=a["probe_selection"],
        lr_mirror=a["lr_mirror"], missing=a["missing"], norm_matched=a["norm_matched"],
        donors=a["donors"])
    expr.to_csv(out)
    print(f"[write] {out}  shape={expr.shape}")
    return out


def main() -> None:
    cfg = load_config()
    ap = argparse.ArgumentParser()
    ap.add_argument("--parc", required=True, choices=["schaefer100", "yan600"])
    ap.add_argument("--pheno", default=cfg["gwas"]["primary"])
    ap.add_argument("--n-spins", type=int, default=cfg["enrichment"]["n_spins"])
    args = ap.parse_args()
    seed = cfg["seed"]

    glab, info = global_labels(args.parc)
    expr_path = build_and_run_abagen(args.parc, glab, info, cfg)

    ids, cents, hemi, target = parc_centroids_and_mask(glab, info)
    expr = pd.read_csv(expr_path, index_col=0).loc[ids]
    geneZ = pd.read_csv(DERIVED / f"{args.pheno}_geneZ.csv")
    top_genes = (DERIVED / f"{args.pheno}_genes_top100.txt").read_text().split()

    scores = {"continuous": EL.continuous_score(expr, geneZ, cfg["enrichment"]["corr"]),
              "top100": EL.topn_score(expr, top_genes)}

    spin_path = RESULTS / "derived" / f"spins_{args.parc}_{args.n_spins}_s{seed}.npy"
    spin_path.parent.mkdir(parents=True, exist_ok=True)
    if spin_path.exists():
        spins = np.load(spin_path)
    else:
        spins = EL.gen_spins(cents, hemi, args.n_spins, seed)
        np.save(spin_path, spins)

    rows = []
    for sname, svec in scores.items():
        res = EL.spin_pvalue(np.asarray(svec), target, spins, cfg["enrichment"]["alternative"])
        res.update({"parcellation": args.parc, "n_parcels": len(ids),
                    "n_target": int(target.sum()), "phenotype": args.pheno, "score": sname})
        rows.append(res)
        print(f"  {args.parc} {sname:10s}: n_target={int(target.sum())}/{len(ids)} "
              f"diff={res['diff']:+.4f} z={res['z_vs_null']:+.2f} p_spin={res['p_spin']:.4f}")

    df = pd.DataFrame(rows)[["parcellation", "n_parcels", "phenotype", "score",
                             "n_target", "diff", "z_vs_null", "p_spin"]]
    out = RESULTS / "tables" / "enrichment_parcellations.csv"
    if out.exists():
        prev = pd.read_csv(out)
        prev = prev[~((prev.parcellation == args.parc) & (prev.phenotype == args.pheno))]
        df = pd.concat([prev, df], ignore_index=True)
    df.to_csv(out, index=False)
    print(f"[write] {out}")

    log_run("06_parcellation_sensitivity", {
        "parcellation": args.parc, "n_parcels": int(len(ids)),
        "phenotype": args.pheno, "overlap_thr": OVERLAP_THR,
        "results": rows})


if __name__ == "__main__":
    main()
