"""Fase 1 (post-process) — annotate MAGMA gene-based output & export gene sets.

Reads results/tables/magma/<pheno>.genes.out[.txt] (the Windows MAGMA build
appends .txt), maps Entrez gene IDs -> HGNC symbols using the NCBI37.3.gene.loc
symbol column (offline, reproducible), applies multiple-testing correction, and
writes:
  results/tables/<pheno>_genes_ranked.csv   full table (symbol,Z,P,P_bonf,P_fdr)
  data/derived/<pheno>_geneZ.csv            symbol,ZSTAT,P for rank/continuous use
  data/derived/<pheno>_genes_<set>.txt      candidate gene sets (see SETS)
  data/derived/alpha_genes_magma.txt        alias of the config-selected set for
                                            the primary phenotype (used by Fase 3)

Run:  python scripts/01_postprocess_magma.py --pheno alphaCz
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from statsmodels.stats.multitest import multipletests

from utils import ROOT, load_config, log_run

DERIVED = ROOT / "data" / "derived"
GENE_LOC = ROOT / "data" / "reference" / "NCBI37.3.gene.loc"

# Candidate gene-set definitions (name -> selector on the ranked frame).
SETS = {
    "bonferroni": lambda d: d["P_bonf"] < 0.05,
    "fdr05":      lambda d: d["P_fdr"] < 0.05,
    "p001":       lambda d: d["P"] < 1e-3,
    "p01":        lambda d: d["P"] < 1e-2,
    "top100":     lambda d: d["P"].rank(method="first") <= 100,
    "top200":     lambda d: d["P"].rank(method="first") <= 200,
}


def resolve_genes_out(pheno: str) -> Path:
    base = ROOT / "results" / "tables" / "magma" / f"{pheno}.genes.out"
    for cand in (base, base.with_suffix(".out.txt"), Path(str(base) + ".txt")):
        if cand.exists():
            return cand
    raise SystemExit(
        f"MAGMA output not found for '{pheno}'. Run "
        f"`python scripts/01_gene_based_magma.py --pheno {pheno}` first."
    )


def entrez_to_symbol() -> dict:
    gl = pd.read_csv(GENE_LOC, sep=r"\s+", header=None,
                     names=["entrez", "chr", "start", "stop", "strand", "symbol"],
                     dtype=str)
    return dict(zip(gl["entrez"], gl["symbol"]))


def main() -> None:
    cfg = load_config()
    ap = argparse.ArgumentParser()
    ap.add_argument("--pheno", default=cfg["gwas"]["primary"])
    args = ap.parse_args()

    magma_out = resolve_genes_out(args.pheno)
    df = pd.read_csv(magma_out, sep=r"\s+")
    df["symbol"] = df["GENE"].astype(str).map(entrez_to_symbol()).fillna(df["GENE"].astype(str))
    df["P_bonf"] = multipletests(df["P"], method="bonferroni")[1]
    df["P_fdr"] = multipletests(df["P"], method="fdr_bh")[1]
    df = df.sort_values("P").reset_index(drop=True)

    ranked = ROOT / "results" / "tables" / f"{args.pheno}_genes_ranked.csv"
    df.to_csv(ranked, index=False)

    # full gene-level Z/P (for rank-based / continuous enrichment, no threshold)
    geneZ = DERIVED / f"{args.pheno}_geneZ.csv"
    df[["symbol", "ZSTAT", "P"]].to_csv(geneZ, index=False)

    # candidate gene sets
    counts = {}
    for name, sel in SETS.items():
        genes = df.loc[sel(df), "symbol"].dropna().unique().tolist()
        (DERIVED / f"{args.pheno}_genes_{name}.txt").write_text(
            "\n".join(genes) + ("\n" if genes else ""))
        counts[name] = len(genes)

    # alias the config-selected set for the primary phenotype (Fase 3 default input)
    chosen = cfg["gwas"].get("gene_set", "top100")
    if args.pheno == cfg["gwas"]["primary"]:
        src = DERIVED / f"{args.pheno}_genes_{chosen}.txt"
        (DERIVED / "alpha_genes_magma.txt").write_text(src.read_text())

    print(f"[write] {ranked}  ({len(df)} genes)")
    print(f"[write] {geneZ}")
    print("gene-set sizes:", counts)
    print(f"top 5: {df.loc[:4, ['symbol','CHR','ZSTAT','P','P_fdr']].to_dict('records')}")
    print(f"primary alias -> '{chosen}' set ({counts.get(chosen)} genes)")

    log_run("01_postprocess_magma", {
        "phenotype": args.pheno, "n_genes": int(len(df)),
        "gene_set_sizes": counts, "selected_set": chosen,
        "magma_out": str(magma_out.relative_to(ROOT))})


if __name__ == "__main__":
    main()
