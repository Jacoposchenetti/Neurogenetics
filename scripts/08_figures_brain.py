"""Brain-surface figures (fsaverage5 inflated): Fig 2 (atlas ROIs) + Fig 4 (score).

Fig 2  Alpha-source ROIs on the inflated cortex for Glasser-360, Schaefer-100 and
       Yan-600. All atlas parcel boundaries drawn as thin lines; only the target
       parcels filled (Glasser by sector; Schaefer/Yan in the alpha accent).
Fig 4  The continuous alpha transcriptomic score painted on the surface with the
       41 Glasser alpha-source ROIs outlined.

Run:  python scripts/08_figures_brain.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import nibabel as nib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, TwoSlopeNorm, LinearSegmentedColormap
from matplotlib.patches import Patch
from nilearn import datasets, plotting, surface

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import ROOT, load_config, log_run
import enrichment_lib as EL
sys.path.insert(0, str(ROOT / "paper" / "figures"))
import _style as S

REFERENCE = ROOT / "data" / "reference"
FIGDIR = ROOT / "paper" / "figures"
DERIVED = ROOT / "data" / "derived"

ROI_BY_SECTOR = {
    "occipital": ["V1", "V2", "ProS", "V3", "V4", "V6", "V6A", "V7", "IPS1", "V3A",
                  "V3B", "V3CD", "IP0", "PGp", "LO1", "LO2"],
    "parietal": ["1", "2", "3a", "3b", "4", "6mp", "6d"],
    "frontal": ["8BL", "9p", "9m", "9a", "8Ad", "9-46d", "8BM", "8Av", "46", "8C",
                "p9-46v", "a32pr", "d32", "a9-46v", "10d", "p10p", "p47r", "IFSa"],
}
VIEWS = [("left", "lateral"), ("left", "medial"), ("right", "lateral"), ("right", "medial")]


def fs5():
    return datasets.fetch_surf_fsaverage("fsaverage5")


def mesh_neighbors(mesh_file):
    coords, faces = surface.load_surf_mesh(mesh_file)
    return faces


def boundary_mask(labels, faces):
    """Vertices lying on a parcel boundary (a face-neighbour has a different label)."""
    b = np.zeros(labels.shape[0], dtype=bool)
    for a, c in ((0, 1), (1, 2), (0, 2)):
        e0, e1 = faces[:, a], faces[:, c]
        diff = labels[e0] != labels[e1]
        b[e0[diff]] = True
        b[e1[diff]] = True
    return b


def glasser_vertex(hemi):
    lab = np.asarray(nib.load(REFERENCE / f"glasser_fsaverage5_{'lh' if hemi=='left' else 'rh'}.label.gii").agg_data()).astype(int)
    return lab


def alpha_source_ids():
    info = pd.read_csv(REFERENCE / "glasser_360_info.csv")
    flat = [l for v in ROI_BY_SECTOR.values() for l in v]
    return info[info["label"].isin(flat)].set_index("id")["label"].to_dict()


def sector_of(label):
    for sec, labs in ROI_BY_SECTOR.items():
        if label in labs:
            return sec
    return None


SECTOR_CODE = {"occipital": 1, "parietal": 2, "frontal": 3}


def glasser_sector_vertex(hemi):
    """Per-vertex sector code from the Glasser alpha-source ROIs (0=none,1/2/3)."""
    lab = glasser_vertex(hemi)
    id2lab = alpha_source_ids()
    sec = np.zeros(lab.shape[0], dtype=int)
    for vid, label in id2lab.items():
        sec[lab == vid] = SECTOR_CODE[sector_of(label)]
    return sec


# ------- per-vertex code map: 0 cortex, 1 occ, 2 par, 3 front, 4 boundary, 5 medial -------
CODE_COLORS = [S.CORTEX_GREY, S.SECTOR["occipital"], S.SECTOR["parietal"],
               S.SECTOR["frontal"], S.BOUNDARY, "#ffffff"]
CODE_CMAP = ListedColormap(CODE_COLORS)


def atlas_target_codemap(atlas, hemi, faces):
    """Return per-vertex code map for Fig 2."""
    if atlas == "glasser":
        lab = glasser_vertex(hemi)
        id2lab = alpha_source_ids()
        sec_code = {1: "occipital", 2: "parietal", 3: "frontal"}
        code = np.zeros(lab.shape[0], dtype=int)
        for vid in np.unique(lab[lab != 0]):
            if vid in id2lab:
                sec = sector_of(id2lab[vid])
                c = {"occipital": 1, "parietal": 2, "frontal": 3}[sec]
                code[lab == vid] = c
        full_labels = lab
    else:
        if atlas == "yan600":
            annot = REFERENCE / f"yan600_{'lh' if hemi=='left' else 'rh'}.annot"
        else:  # schaefer100
            from netneurotools.datasets import fetch_schaefer2018
            s = fetch_schaefer2018(version="fsaverage5")["100Parcels7Networks"]
            annot = s[0] if hemi == "left" else s[1]
        lab = nib.freesurfer.read_annot(str(annot))[0].astype(int)
        # per-vertex Glasser sector code (1 occ, 2 par, 3 front; 0 = not a source)
        sec_v = glasser_sector_vertex(hemi)
        code = np.zeros(lab.shape[0], dtype=int)
        for pid in np.unique(lab[lab != 0]):
            v = lab == pid
            srcs = sec_v[v]
            if (srcs > 0).mean() >= 0.5:            # target parcel
                code[v] = np.bincount(srcs[srcs > 0]).argmax()  # majority sector
        full_labels = lab
    # boundaries of ALL parcels
    b = boundary_mask(full_labels, faces)
    code[b & (full_labels != 0)] = 4
    code[full_labels == 0] = 5  # medial wall
    return code.astype(float)


def fig2_atlases():
    fs = fs5()
    atlases = [("glasser", "Glasser HCP-MMP1.0 (360)"),
               ("schaefer100", "Schaefer (100)"),
               ("yan600", "Yan homotopic (600)")]
    fig, axes = plt.subplots(len(atlases), 4, figsize=(13, 8.2),
                             subplot_kw={"projection": "3d"})
    for r, (atlas, title) in enumerate(atlases):
        # for schaefer/yan the accent colour replaces the sector palette
        for c, (hemi, view) in enumerate(VIEWS):
            mesh = fs[f"infl_{'left' if hemi=='left' else 'right'}"]
            bg = fs[f"sulc_{'left' if hemi=='left' else 'right'}"]
            code = atlas_target_codemap(atlas, hemi, mesh_neighbors(mesh))
            cmap = CODE_CMAP   # sector colours (occ/par/front) for every atlas
            plotting.plot_surf_roi(mesh, roi_map=code, hemi=hemi, view=view,
                                   bg_map=bg, darkness=0.5, cmap=cmap,
                                   vmin=0, vmax=5, axes=axes[r, c], figure=fig,
                                   colorbar=False)
            if r == 0:
                axes[r, c].set_title(f"{hemi[:1].upper()}H {view}", fontsize=10)
        axes[r, 0].text2D(-0.08, 0.5, title, transform=axes[r, 0].transAxes,
                          rotation=90, va="center", ha="center", fontsize=11, color=S.INK)
    handles = [Patch(fc=S.SECTOR["occipital"], label="Occipital ROIs"),
               Patch(fc=S.SECTOR["parietal"], label="Parietal ROIs"),
               Patch(fc=S.SECTOR["frontal"], label="Frontal ROIs"),
               Patch(fc=S.CORTEX_GREY, label="Non-target cortex")]
    fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=9, bbox_to_anchor=(0.5, 0.005))
    fig.suptitle("Alpha-source ROIs (Tabarelli et al. 2022) on the inflated cortex, across parcellations",
                 fontsize=13)
    fig.tight_layout(rect=[0.02, 0.05, 1, 0.97])
    out = FIGDIR / "fig2_atlas_rois.png"
    fig.savefig(out); plt.close(fig)
    print(f"[write] {out}")


def _autocrop(img, pad=8):
    """Trim near-white margins from an RGB(A) image array."""
    rgb = img[..., :3] if img.shape[-1] == 4 else img
    mask = (rgb < 0.985).any(axis=-1)
    if not mask.any():
        return img
    rows, cols = np.where(mask.any(1))[0], np.where(mask.any(0))[0]
    r0, r1 = max(rows.min() - pad, 0), min(rows.max() + pad, img.shape[0])
    c0, c1 = max(cols.min() - pad, 0), min(cols.max() + pad, img.shape[1])
    return img[r0:r1, c0:c1]


def fig4_scoremap():
    """Vivid continuous score on the surface, rendered with the plotly engine
    (matplotlib's 3D shading mutes continuous colours; plotly does not). ROI
    boundaries are baked in as background (grey) lines by NaN-masking boundary
    vertices. Panels are exported to PNG (kaleido) and assembled in matplotlib."""
    import matplotlib.image as mpimg
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import Normalize
    cfg = load_config()
    ids, cents, hemi = EL.parcel_centroids()
    expr = pd.read_csv(DERIVED / "ahba_expression_glasser.csv", index_col=0).loc[ids]
    geneZ = pd.read_csv(DERIVED / "alphaCz_geneZ.csv")
    raw = pd.Series(EL.continuous_score(expr, geneZ, cfg["enrichment"]["corr"]), index=ids)
    score = (raw - raw.mean()) / raw.std()
    id2lab = alpha_source_ids()
    fs = fs5()
    vmax = 2.5
    tmpdir = Path(__file__).resolve().parents[1]
    import tempfile
    tmp = Path(tempfile.mkdtemp())
    panels = []
    for hemi_, view in VIEWS:
        mesh = fs[f"infl_{'left' if hemi_=='left' else 'right'}"]
        bg = fs[f"sulc_{'left' if hemi_=='left' else 'right'}"]
        lab = glasser_vertex(hemi_)
        vmap = np.full(lab.shape[0], np.nan)
        for vid in np.unique(lab[lab != 0]):
            if vid in score.index:
                vmap[lab == vid] = score.loc[vid]
        # thin ROI outline: only the OUTER boundary (target vertices adjacent to
        # non-target), not internal target-target seams -> a single thin ring.
        is_t = np.isin(lab, list(id2lab))
        outline = boundary_mask(is_t.astype(int), mesh_neighbors(mesh)) & is_t
        vmap[outline] = np.nan
        f = plotting.plot_surf(mesh, surf_map=vmap, hemi=hemi_, view=view, bg_map=bg,
                               bg_on_data=True, darkness=0.3, cmap="RdBu_r",
                               vmin=-vmax, vmax=vmax, engine="plotly", symmetric_cmap=True)
        # zoom the camera out so the whole brain fits (nilearn's default eye clips it)
        eye = f.figure.layout.scene.camera.eye
        f.figure.update_layout(scene_camera=dict(eye=dict(x=eye.x * 1.6, y=eye.y * 1.6,
                                                          z=eye.z * 1.6)),
                               margin=dict(l=0, r=0, t=0, b=0))
        p = tmp / f"{hemi_}_{view}.png"
        f.figure.write_image(str(p), width=760, height=760, scale=3)  # square hi-res
        panels.append((f"{hemi_[:1].upper()}H {view}", mpimg.imread(str(p))))

    fig, axes = plt.subplots(2, 2, figsize=(9, 8.4))   # square 2x2: more room per brain
    for ax, (title, img) in zip(axes.ravel(), panels):
        ax.imshow(_autocrop(img)); ax.set_title(title, fontsize=11); ax.axis("off")
    sm = ScalarMappable(norm=Normalize(-vmax, vmax), cmap="RdBu_r"); sm.set_array([])
    cax = fig.add_axes([0.36, 0.055, 0.28, 0.022])
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal", ticks=[-vmax, 0, vmax])
    cb.set_label("alpha transcriptomic score (z across regions)", fontsize=9)
    cb.ax.set_xticklabels(["low", "0", "high"])
    fig.suptitle("Continuous alpha transcriptomic score (alphaCz) on cortex\ndark outline = 41 alpha-source ROIs",
                 fontsize=12)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.90, bottom=0.09, wspace=0.02, hspace=0.08)
    out = FIGDIR / "fig4_score_map.png"
    fig.savefig(out, dpi=200); plt.close(fig)
    print(f"[write] {out}")


def main():
    S.apply_rc()
    fig2_atlases()
    fig4_scoremap()
    log_run("08_figures_brain", {"figures": ["fig2_atlas_rois.png", "fig4_score_map.png"]})


if __name__ == "__main__":
    main()
