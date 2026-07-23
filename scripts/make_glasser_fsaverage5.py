"""Resample the Glasser HCP-MMP1.0 atlas from fsLR32k to fsaverage5 (pip-only).

Why: Tabarelli et al. 2022 defines alpha sources as Glasser parcels. Glasser is
native to fsLR32k, but abagen maps AHBA samples to fsaverage5, so the atlas must
be delivered as a tuple of fsaverage5 GIFTI label images.

How (no Connectome Workbench needed): fsLR and fsaverage use DIFFERENT spherical
registrations, so a naive sphere nearest-neighbor is wrong. We instead use the
canonical HCP `resample_fsaverage` sphere pair, in which BOTH surfaces live in the
common *fsaverage* spherical frame:
  - fs_LR-deformed_to-fsaverage.{L,R}.sphere.32k_fs_LR  : the 32k fsLR vertices,
    warped into the fsaverage frame (vertex order matches the mmpall labels).
  - fsaverage5_std_sphere.{L,R}.10k_fsavg                : the fsaverage5 (10k)
    vertices in the same frame (vertex order matches FreeSurfer fsaverage5, which
    is what abagen expects).
For each fsaverage5 vertex we take the label of the nearest fsLR vertex (KDTree in
the shared frame). This is exactly the correspondence Workbench uses; nearest-
neighbor is the standard choice for discrete label resampling.

Outputs (into data/reference/):
  glasser_fsaverage5_lh.label.gii, glasser_fsaverage5_rh.label.gii
  glasser_360_info.csv   (abagen atlas_info: id,label,hemisphere,structure)

Run:  python scripts/make_glasser_fsaverage5.py
Trusted sources (OSF, via netneurotools): mmpall labels + HCP standard meshes.
"""
from __future__ import annotations

import os
import re
from pathlib import Path

import numpy as np
import nibabel as nib
import pandas as pd
from scipy.spatial import cKDTree

from netneurotools.datasets import fetch_mmpall, fetch_hcp_standards

from utils import ROOT, log_run

REFERENCE = ROOT / "data" / "reference"
RESAMPLE_SUBDIR = "resample_fsaverage"


def _pointset(gii_path: str) -> np.ndarray:
    """Return the (n_vertices, 3) coordinate array of a .surf.gii."""
    return np.asarray(nib.load(gii_path).agg_data("NIFTI_INTENT_POINTSET"))


def _strip_label(raw: str) -> tuple[str, str]:
    """'L_8C_ROI' -> ('8C', 'L'); '???' -> ('', '')."""
    m = re.match(r"^([LR])_(.+)_ROI$", raw)
    if not m:
        return "", ""
    hemi, name = m.group(1), m.group(2)
    return name, hemi


def resample_hemi(label_gii: str, fslr_sphere: str, fsavg5_sphere: str) -> np.ndarray:
    labels = np.asarray(nib.load(label_gii).agg_data()).astype(np.int32)  # (32492,)
    fslr_xyz = _pointset(fslr_sphere)                                     # (32492, 3)
    fsavg5_xyz = _pointset(fsavg5_sphere)                                 # (10242, 3)
    if labels.shape[0] != fslr_xyz.shape[0]:
        raise ValueError(
            f"vertex mismatch: labels {labels.shape[0]} vs fsLR sphere {fslr_xyz.shape[0]}"
        )
    # nearest fsLR vertex (in fsaverage frame) for each fsaverage5 vertex
    _, idx = cKDTree(fslr_xyz).query(fsavg5_xyz, k=1)
    return labels[idx]  # (10242,)


def build_atlas_info(label_tables) -> pd.DataFrame:
    """id,label,hemisphere,structure for all 360 parcels (from a mmpall labeltable)."""
    rows = []
    for lab in label_tables:
        if lab.key == 0:
            continue
        name, hemi = _strip_label(lab.label)
        if not name:
            continue
        rows.append({"id": int(lab.key), "label": name,
                     "hemisphere": hemi, "structure": "cortex"})
    df = pd.DataFrame(rows).sort_values("id").reset_index(drop=True)
    return df


def write_label_gifti(labels: np.ndarray, info: pd.DataFrame, out_path: Path) -> None:
    """Write an fsaverage5 integer-label GIFTI with an embedded label table."""
    lt = nib.gifti.GiftiLabelTable()
    bg = nib.gifti.GiftiLabel(key=0, red=0, green=0, blue=0, alpha=0)
    bg.label = "???"
    lt.labels.append(bg)
    for _, r in info.iterrows():
        gl = nib.gifti.GiftiLabel(key=int(r["id"]), red=0.5, green=0.5, blue=0.5, alpha=1)
        gl.label = f'{r["hemisphere"]}_{r["label"]}'
        lt.labels.append(gl)
    darr = nib.gifti.GiftiDataArray(
        labels.astype(np.int32), intent="NIFTI_INTENT_LABEL", datatype="NIFTI_TYPE_INT32"
    )
    gii = nib.gifti.GiftiImage(labeltable=lt, darrays=[darr])
    gii.to_filename(str(out_path))


def main() -> None:
    REFERENCE.mkdir(parents=True, exist_ok=True)
    mmp = fetch_mmpall()                       # .L / .R fsLR32k label giftis
    std = Path(fetch_hcp_standards()) / RESAMPLE_SUBDIR

    info = build_atlas_info(nib.load(mmp.L).labeltable.labels)
    info_path = REFERENCE / "glasser_360_info.csv"
    info.to_csv(info_path, index=False)

    hemis = {
        "lh": dict(label=mmp.L,
                   fslr=std / "fs_LR-deformed_to-fsaverage.L.sphere.32k_fs_LR.surf.gii",
                   fsavg5=std / "fsaverage5_std_sphere.L.10k_fsavg_L.surf.gii"),
        "rh": dict(label=mmp.R,
                   fslr=std / "fs_LR-deformed_to-fsaverage.R.sphere.32k_fs_LR.surf.gii",
                   fsavg5=std / "fsaverage5_std_sphere.R.10k_fsavg_R.surf.gii"),
    }

    summary = {}
    for hemi, f in hemis.items():
        lab5 = resample_hemi(str(f["label"]), str(f["fslr"]), str(f["fsavg5"]))
        out = REFERENCE / f"glasser_fsaverage5_{hemi}.label.gii"
        write_label_gifti(lab5, info, out)
        n_nonzero = int(np.count_nonzero(lab5))
        n_parcels = int(np.unique(lab5[lab5 != 0]).size)
        summary[hemi] = {"n_vertices": int(lab5.size), "n_cortex_vertices": n_nonzero,
                         "n_parcels": n_parcels, "output": str(out.relative_to(ROOT))}
        print(f"[{hemi}] {lab5.size} verts, {n_nonzero} cortical, {n_parcels} parcels -> {out.name}")
        # Glasser has exactly 180 areas per hemisphere; all should survive resampling
        if n_parcels != 180:
            print(f"  [WARN] expected 180 parcels, got {n_parcels}")

    print(f"[write] {info_path}  ({len(info)} parcels)")
    log_run("make_glasser_fsaverage5",
            {"n_parcels_total": int(len(info)), "hemispheres": summary,
             "source": "netneurotools fetch_mmpall + fetch_hcp_standards (OSF)",
             "method": "nearest-neighbor in fsaverage frame via resample_fsaverage sphere pair"})


if __name__ == "__main__":
    main()
