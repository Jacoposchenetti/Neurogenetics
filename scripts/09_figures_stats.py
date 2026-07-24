"""Statistical figures: Fig 1 (pipeline), Fig 3 (genetics), Fig 5 (enrichment),
Fig 6 (validation). Standard matplotlib; palette from paper/figures/_style.py.

Run:  python scripts/09_figures_stats.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import ROOT, load_config, log_run
import enrichment_lib as EL
sys.path.insert(0, str(ROOT / "paper" / "figures"))
import _style as S

DERIVED = ROOT / "data" / "derived"
TABLES = ROOT / "results" / "tables"
FIGDIR = ROOT / "paper" / "figures"


def gene_table(pheno: str = "alphaCz") -> pd.DataFrame:
    """Gene-level results with symbols, from the versioned ranked table.

    Uses results/tables/<pheno>_genes_ranked.csv, which 01_postprocess_magma.py
    already wrote with the Entrez->symbol mapping resolved (columns GENE, CHR,
    START, symbol, ZSTAT, P, ...). Reading it here keeps the figures runnable
    from the committed derived results alone: neither the MAGMA .genes.out nor
    data/reference/NCBI37.3.gene.loc (part of the non-redistributed reference
    panel) is needed to rebuild Fig 3.
    """
    ranked = TABLES / f"{pheno}_genes_ranked.csv"
    if not ranked.exists():
        raise SystemExit(
            f"missing {ranked.relative_to(ROOT)} — run "
            f"`python scripts/01_postprocess_magma.py --pheno {pheno}` first."
        )
    return pd.read_csv(ranked)


# --------------------------------- Fig 1 pipeline ---------------------------------
def fig1_pipeline():
    S.apply_rc()
    fig, ax = plt.subplots(figsize=(12, 3.2)); ax.axis("off"); ax.set_xlim(0, 12); ax.set_ylim(0, 3)
    steps = [
        ("ENIGMA-EEG\nGWAS\n(6 phenotypes)", S.OKABE["blue"]),
        ("MAGMA\ngene-based\n(gene-level Z)", S.OKABE["skyblue"]),
        ("AHBA / abagen\nregional expression\n(Glasser 360)", S.OKABE["green"]),
        ("Alpha generators\nTabarelli 2022\n(41 ROIs)", S.OKABE["orange"]),
        ("Enrichment vs\nspatial + gene-set\n+ control nulls", S.OKABE["vermillion"]),
    ]
    w, h, y = 2.0, 1.4, 1.4
    xs = np.linspace(0.4, 9.6, len(steps))
    for i, (txt, col) in enumerate(steps):
        x = xs[i]
        box = FancyBboxPatch((x, y - h / 2), w, h, boxstyle="round,pad=0.05,rounding_size=0.12",
                             fc=col, ec="none", alpha=0.9)
        ax.add_patch(box)
        ax.text(x + w / 2, y, txt, ha="center", va="center", color="white", fontsize=9.5, weight="bold")
        if i < len(steps) - 1:
            ax.add_patch(FancyArrowPatch((x + w, y), (xs[i + 1], y), arrowstyle="-|>",
                                         mutation_scale=16, color=S.INK, lw=1.6))
    ax.text(6, 2.8, "Do EEG alpha-power genes concentrate in the cortical generators of the alpha rhythm?",
            ha="center", va="center", fontsize=11, style="italic", color=S.INK)
    fig.tight_layout(); out = FIGDIR / "fig1_pipeline.png"; fig.savefig(out); plt.close(fig)
    print(f"[write] {out}")


# --------------------------------- Fig 3 genetics ---------------------------------
def fig3_genetics():
    S.apply_rc()
    df = gene_table("alphaCz")
    df["logp"] = -np.log10(df["P"])
    # cumulative genomic position
    df = df[df["CHR"].apply(lambda x: str(x).isdigit())].copy()
    df["CHR"] = df["CHR"].astype(int); df = df.sort_values(["CHR", "START"])
    offs, cum, ticks, ticklab = {}, 0, [], []
    for ch in range(1, 23):
        sub = df[df["CHR"] == ch]
        if sub.empty:
            continue
        offs[ch] = cum
        ticks.append(cum + sub["START"].max() / 2); ticklab.append(str(ch))
        cum += sub["START"].max() + 2e7
    df["pos"] = df.apply(lambda r: offs[r["CHR"]] + r["START"], axis=1)

    fig = plt.figure(figsize=(12, 6.2))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.1, 1], hspace=0.42, wspace=0.28)
    # (A) Manhattan
    axA = fig.add_subplot(gs[0, :])
    for ch in range(1, 23):
        sub = df[df["CHR"] == ch]
        axA.scatter(sub["pos"], sub["logp"], s=7,
                    color=(S.OKABE["blue"] if ch % 2 else "#9aa0a8"), rasterized=True)
    bonf = -np.log10(0.05 / len(df))
    axA.axhline(bonf, ls="--", lw=1, color=S.OKABE["vermillion"], label=f"Bonferroni (p={0.05/len(df):.1e})")
    axA.set_ylim(0, 7.4)   # headroom so the title clears PRKG2 and labels clear the line
    prkg2 = df[df["symbol"] == "PRKG2"].iloc[0]
    axA.annotate("PRKG2", (prkg2["pos"], prkg2["logp"]), fontsize=8.5, ha="center",
                 va="bottom", xytext=(0, 4), textcoords="offset points", color=S.INK, weight="bold")
    clust = df[df["symbol"].isin(["GNL3", "GLT8D1", "NT5DC2", "PBRM1", "SPCS1"])]
    if not clust.empty:
        cx, cy = clust["pos"].mean(), clust["logp"].max()
        # label in the empty band (below Bonferroni, above the data cloud), arrow to cluster
        axA.annotate("chr3p21 cluster\n(GNL3, GLT8D1, NT5DC2, PBRM1)", xy=(cx, cy),
                     xytext=(ticks[5], 5.05), textcoords="data", fontsize=7.5,
                     ha="center", va="center", color=S.INK,
                     arrowprops=dict(arrowstyle="->", lw=0.7, color=S.INK,
                                     connectionstyle="arc3,rad=-0.15"))
    axA.set_xticks(ticks); axA.set_xticklabels(ticklab, fontsize=7)
    axA.set_ylabel("gene-based  −log₁₀(p)"); axA.set_xlabel("chromosome")
    axA.set_title("A  MAGMA gene-based association, central alpha power (alphaCz)", pad=14)
    axA.legend(fontsize=8, loc="upper right"); axA.grid(axis="x", visible=False)
    # (B) top-10 bar
    axB = fig.add_subplot(gs[1, 0])
    top = df.nsmallest(10, "P")[::-1]
    axB.barh(top["symbol"], top["logp"], color=S.OKABE["blue"], height=0.62)
    axB.set_xlabel("−log₁₀(p)"); axB.set_title("B  Top 10 genes (alphaCz)"); axB.grid(axis="y", visible=False)
    # (C) top genes across phenotypes
    axC = fig.add_subplot(gs[1, 1])
    phs = ["alphaCz", "alphaOcc", "peakOcc", "thetaCz", "betaCz", "deltaCz"]
    topsyms = df.nsmallest(8, "P")["symbol"].tolist()
    mat = np.full((len(topsyms), len(phs)), np.nan)
    for j, ph in enumerate(phs):
        g = pd.read_csv(DERIVED / f"{ph}_geneZ.csv").drop_duplicates("symbol").set_index("symbol")
        for i, sym in enumerate(topsyms):
            if sym in g.index:
                mat[i, j] = -np.log10(g.loc[sym, "P"])
    im = axC.imshow(mat, aspect="auto", cmap="cividis", vmin=0)
    axC.set_xticks(range(len(phs))); axC.set_xticklabels(phs, rotation=45, ha="right", fontsize=8)
    axC.set_yticks(range(len(topsyms))); axC.set_yticklabels(topsyms, fontsize=8)
    axC.set_title("C  Top alphaCz genes across phenotypes"); axC.grid(False)
    fig.colorbar(im, ax=axC, fraction=0.046, pad=0.04, label="−log₁₀(p)")
    out = FIGDIR / "fig3_genetics.png"; fig.savefig(out, bbox_inches="tight"); plt.close(fig)
    print(f"[write] {out}")


# ------------------------- Fig 5 enrichment / specificity -------------------------
def fig5_enrichment():
    S.apply_rc()
    cfg = load_config()
    spins, ids = EL.load_or_make_spins(cfg["enrichment"]["n_spins"], cfg["seed"])
    masks = EL.masks_for(ids)
    scores = EL.score_maps("alphaCz", ids, cfg["enrichment"]["corr"])
    allp = pd.read_csv(TABLES / "enrichment_all_phenotypes.csv")
    parc = pd.read_csv(TABLES / "enrichment_parcellations.csv")

    fig = plt.figure(figsize=(12, 7.4))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1], hspace=0.45, wspace=0.28)
    # (A) primary spin nulls
    axA = fig.add_subplot(gs[0, 0])
    for sname, col in [("continuous", S.OKABE["vermillion"]), ("top100", S.OKABE["blue"])]:
        sv = np.asarray(scores[sname]); m = masks["full"]
        obs = sv[m].mean()
        null = np.array([sv[spins[:, s]][m].mean() for s in range(spins.shape[1])])
        p = (1 + np.sum(null >= obs)) / (len(null) + 1)
        axA.hist(null, bins=50, color=col, alpha=0.35, density=True)
        axA.axvline(obs, color=col, lw=2, label=f"{sname}: p={p:.3f}")
    axA.set_title("A  Primary enrichment vs spin null (alphaCz, 41 ROIs)")
    axA.set_xlabel("mean alpha score in target"); axA.set_ylabel("density"); axA.legend(fontsize=8)
    # (B) specificity bars
    axB = fig.add_subplot(gs[0, 1])
    sub = allp[allp.roi_set == "full"].copy()
    order = ["alphaCz", "alphaOcc", "peakOcc", "thetaCz", "betaCz", "deltaCz"]
    x = np.arange(len(order)); wd = 0.38
    for k, (sc, off) in enumerate([("continuous", -wd / 2), ("top100", wd / 2)]):
        d = sub[sub.score == sc].set_index("phenotype").loc[order]
        cols = [S.ALPHA_ACCENT if b != "control" else S.CONTROL_GREY for b in d.band]
        axB.bar(x + off, -np.log10(d.p_spin), wd, color=cols, alpha=0.85 if sc == "continuous" else 0.55,
                edgecolor="white", label=sc)
    axB.axhline(-np.log10(0.05), ls="--", lw=1, color=S.INK)
    axB.set_xticks(x); axB.set_xticklabels(order, rotation=45, ha="right", fontsize=8)
    axB.set_ylabel("−log₁₀(p_spin)")
    axB.set_title("B  Specificity: alpha (colour) vs controls (grey); nothing survives FDR")
    axB.legend(fontsize=8, title="score", title_fontsize=8); axB.grid(axis="x", visible=False)
    # (C) parcellation sensitivity
    axC = fig.add_subplot(gs[1, 0])
    pmap = {"schaefer100": "Schaefer-100", "glasser360": "Glasser-360", "yan600": "Yan-600"}
    # inject glasser row from primary table
    gl = pd.DataFrame({"parcellation": ["glasser360", "glasser360"], "score": ["continuous", "top100"],
                       "p_spin": [0.021898, 0.030097], "z_vs_null": [1.975, 1.774]})
    pc = pd.concat([parc[["parcellation", "score", "p_spin", "z_vs_null"]], gl], ignore_index=True)
    porder = ["schaefer100", "glasser360", "yan600"]
    x = np.arange(len(porder))
    for sc, off, col in [("continuous", -wd/2, S.OKABE["vermillion"]), ("top100", wd/2, S.OKABE["blue"])]:
        d = pc[pc.score == sc].set_index("parcellation").loc[porder]
        axC.bar(x + off, -np.log10(d.p_spin), wd, color=col, alpha=0.8, edgecolor="white", label=sc)
    axC.axhline(-np.log10(0.05), ls="--", lw=1, color=S.INK)
    axC.set_xticks(x); axC.set_xticklabels([pmap[p] for p in porder], fontsize=9)
    axC.set_ylabel("−log₁₀(p_spin)"); axC.set_title("C  Parcellation sensitivity (alphaCz)")
    axC.legend(fontsize=8); axC.grid(axis="x", visible=False)
    # (D) effect sizes (z) heat-ish table
    axD = fig.add_subplot(gs[1, 1])
    piv = sub.pivot_table(index="phenotype", columns="score", values="z_vs_null").loc[order]
    im = axD.imshow(piv.values, aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5)
    axD.set_xticks(range(2)); axD.set_xticklabels(piv.columns, fontsize=9)
    axD.set_yticks(range(len(order))); axD.set_yticklabels(order, fontsize=8)
    for i in range(len(order)):
        for j in range(2):
            axD.text(j, i, f"{piv.values[i,j]:.2f}", ha="center", va="center", fontsize=8,
                     color="white" if abs(piv.values[i, j]) > 1.5 else S.INK)
    axD.set_title("D  Enrichment effect size (z vs null)"); axD.grid(False)
    fig.colorbar(im, ax=axD, fraction=0.046, pad=0.04, label="z")
    out = FIGDIR / "fig5_enrichment.png"; fig.savefig(out, bbox_inches="tight"); plt.close(fig)
    print(f"[write] {out}")


# ------------------------------ Fig 6 validation ------------------------------
def fig6_validation():
    S.apply_rc()
    cfg = load_config()
    spins, ids = EL.load_or_make_spins(cfg["enrichment"]["n_spins"], cfg["seed"])
    masks = EL.masks_for(ids)
    expr = pd.read_csv(DERIVED / "ahba_expression_glasser.csv", index_col=0).loc[ids]
    rng = np.random.default_rng(cfg["seed"])

    fig, axes = plt.subplots(1, 3, figsize=(12, 3.6))
    # (A) positive control: PC1 gradient spin null
    E = expr.to_numpy(); Ec = E - E.mean(0, keepdims=True)
    U, Sv, Vt = np.linalg.svd(Ec, full_matrices=False)
    pc1_region = U[:, 0] * Sv[0]; pc1_load = Vt[0, :]
    if np.corrcoef(pc1_region, E @ pc1_load)[0, 1] < 0:
        pc1_region, pc1_load = -pc1_region, -pc1_load
    tgt = pc1_region >= np.quantile(pc1_region, 2/3)
    gz = pd.DataFrame({"symbol": expr.columns, "ZSTAT": pc1_load, "P": np.nan})
    sc = np.asarray(EL.continuous_score(expr, gz, cfg["enrichment"]["corr"]))
    obs = sc[tgt].mean() - sc[~tgt].mean()
    null = np.array([sc[spins[:, s]][tgt].mean() - sc[spins[:, s]][~tgt].mean() for s in range(spins.shape[1])])
    p = (1 + np.sum(null >= obs)) / (len(null) + 1)
    axes[0].hist(null, bins=50, color=S.CONTROL_GREY, density=True)
    axes[0].axvline(obs, color=S.OKABE["green"], lw=2, label=f"observed\np={p:.4f}")
    axes[0].set_title("A  Positive control (PC1 gradient)"); axes[0].set_xlabel("target−rest difference")
    axes[0].legend(fontsize=8)
    # (B) gene-set null
    top = [g for g in (DERIVED / "alphaCz_genes_top100.txt").read_text().split() if g in expr.columns]
    def diff_for(gl):
        z = EL.topn_score(expr, list(gl)); return z[masks["full"]].mean() - z[~masks["full"]].mean()
    obs2 = diff_for(top); genes = expr.columns.to_numpy(); k = len(top)
    null2 = np.array([diff_for(genes[rng.choice(len(genes), k, replace=False)]) for _ in range(3000)])
    p2 = (1 + np.sum(null2 >= obs2)) / (len(null2) + 1)
    axes[1].hist(null2, bins=50, color=S.CONTROL_GREY, density=True)
    axes[1].axvline(obs2, color=S.ALPHA_ACCENT, lw=2, label=f"alpha top-100\np={p2:.2f}")
    axes[1].set_title("B  Co-expression gene-set null"); axes[1].set_xlabel("target−rest difference")
    axes[1].legend(fontsize=8)
    # (C) spin vs brainsmash (p-values bar)
    sv = pd.read_csv(TABLES / "validation_spin_validation.csv").iloc[0]
    axes[2].bar(["spin\n(custom)", "brainsmash"], [sv["p_spin_custom"], sv["p_brainsmash"]],
                color=[S.OKABE["blue"], S.OKABE["orange"]], width=0.55)
    axes[2].axhline(0.05, ls="--", lw=1, color=S.INK)
    axes[2].set_ylabel("p-value"); axes[2].set_title("C  Spin validated vs brainsmash")
    for i, v in enumerate([sv["p_spin_custom"], sv["p_brainsmash"]]):
        axes[2].text(i, v + 0.001, f"{v:.3f}", ha="center", fontsize=9)
    fig.tight_layout(); out = FIGDIR / "fig6_validation.png"; fig.savefig(out); plt.close(fig)
    print(f"[write] {out}")


def main():
    fig1_pipeline()
    fig3_genetics()
    fig5_enrichment()
    fig6_validation()
    log_run("09_figures_stats", {"figures": ["fig1_pipeline.png", "fig3_genetics.png",
                                             "fig5_enrichment.png", "fig6_validation.png"]})


if __name__ == "__main__":
    main()
