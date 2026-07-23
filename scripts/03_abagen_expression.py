"""Fase 3 — Regional gene expression from the Allen Human Brain Atlas via abagen.

Produces a region x gene expression matrix in the chosen parcellation and, if a
gene-based hit list from Fase 1 exists, a composite z-scored expression score per
region for the alpha gene set.

Atlas options (--atlas):
  dk       Desikan-Killiany — bundled with abagen, NO external download. Used to
           VALIDATE the pipeline end-to-end and to trigger the AHBA download
           (~4 GB, cached under abagen-data/). NOT the analysis atlas.
  glasser  Glasser HCP-MMP1.0 — the analysis parcellation (matches Tabarelli
           2022). Requires a volumetric MNI152 NIfTI + label table; see
           get_glasser_atlas(). Sourcing is deliberately explicit (no random
           NIfTI from the web).

Run:  python scripts/03_abagen_expression.py --atlas dk
      python scripts/03_abagen_expression.py --atlas glasser
First run downloads AHBA (~4 GB, cached; subsequent runs are fast).

Requires: abagen, nibabel, nilearn, pandas, numpy.  Python 3.12.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from utils import ROOT, load_config, log_run

DERIVED = ROOT / "data" / "derived"
REFERENCE = ROOT / "data" / "reference"


def get_dk_atlas():
    """Desikan-Killiany, bundled with abagen. Returns (atlas_img, atlas_info).

    Zero external download; used to validate the extraction pipeline and to
    trigger/verify the AHBA microarray download.
    """
    import abagen

    dk = abagen.fetch_desikan_killiany()
    return dk["image"], dk["info"]


def get_glasser_atlas():
    """Glasser HCP-MMP1.0 as an fsaverage5 SURFACE atlas. Returns ((lh, rh), info).

    Glasser is native to fsLR32k; abagen maps AHBA samples to fsaverage5, so we
    deliver the atlas as a tuple of fsaverage5 GIFTI label images. These are built
    (pip-only, no Connectome Workbench) by scripts/make_glasser_fsaverage5.py,
    which resamples fsLR32k -> fsaverage5 using the canonical HCP registration
    sphere pair. `atlas_info` labels (e.g. 'V1', '8C') match the Tabarelli ROIs in
    scripts/02_roi_mapping.py.
    """
    lh = REFERENCE / "glasser_fsaverage5_lh.label.gii"
    rh = REFERENCE / "glasser_fsaverage5_rh.label.gii"
    info = REFERENCE / "glasser_360_info.csv"
    if not (lh.exists() and rh.exists() and info.exists()):
        raise FileNotFoundError(
            "fsaverage5 Glasser atlas not found. Build it first:\n"
            "  python scripts/make_glasser_fsaverage5.py\n"
            f"(expected {lh.name}, {rh.name}, {info.name} in data/reference/)."
        )
    return (str(lh), str(rh)), pd.read_csv(info)


ATLASES = {"dk": get_dk_atlas, "glasser": get_glasser_atlas}


def compute_gene_set_score(expression: pd.DataFrame) -> Path | None:
    """Composite z-scored expression score for the Fase 1 alpha gene set."""
    gene_list = DERIVED / "alpha_genes_magma.txt"
    if not gene_list.exists():
        print(f"[skip] {gene_list} not found — run Fase 1 first for the gene-set score.")
        return None
    genes = [g.strip() for g in gene_list.read_text().splitlines() if g.strip()]
    present = [g for g in genes if g in expression.columns]
    if not present:
        print(f"[warn] none of the {len(genes)} alpha genes are in the expression matrix.")
        return None
    sub = expression[present]
    z = (sub - sub.mean(axis=0)) / sub.std(axis=0, ddof=0)
    score = z.mean(axis=1).rename("alpha_expr_score")
    out = DERIVED / "alpha_expr_score.csv"
    score.to_csv(out)
    print(f"[write] {out}  ({len(present)}/{len(genes)} genes matched)")
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--atlas", choices=list(ATLASES), default="dk",
                    help="dk = validate pipeline (bundled); glasser = analysis atlas")
    args = ap.parse_args()

    cfg = load_config()
    import abagen

    atlas_img, atlas_info = ATLASES[args.atlas]()

    a = cfg["abagen"]
    donors = a.get("donors", "all")
    expression = abagen.get_expression_data(
        atlas_img,
        atlas_info=atlas_info,
        ibf_threshold=a["ibf_threshold"],
        probe_selection=a["probe_selection"],
        lr_mirror=a["lr_mirror"],
        missing=a["missing"],
        norm_matched=a["norm_matched"],
        donors=donors,
    )  # -> DataFrame (region_id x gene_symbol)

    out_expr = DERIVED / f"ahba_expression_{args.atlas}.csv"
    expression.to_csv(out_expr)
    print(f"[write] {out_expr}  shape={expression.shape}")

    score_path = compute_gene_set_score(expression)

    log_run(
        "03_abagen_expression",
        {
            "atlas": args.atlas,
            "abagen_params": a,
            "donors": donors,
            "expression_shape": list(expression.shape),
            "n_genes": int(expression.shape[1]),
            "score_output": str(score_path) if score_path else None,
        },
    )


if __name__ == "__main__":
    main()
