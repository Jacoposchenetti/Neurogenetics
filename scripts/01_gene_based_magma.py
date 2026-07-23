"""Fase 1 — MAGMA gene-based analysis from ENIGMA-EEG summary statistics.

Cross-platform driver (Windows-friendly) that:
  1. builds a SNP-location file (SNP CHR BP) from the chosen phenotype sumstats,
  2. runs MAGMA SNP->gene annotation (window-based),
  3. runs the MAGMA gene-based test (SNP-wise, per-SNP N).

Then run scripts/01_postprocess_magma.py to rank/correct and export the gene list.

Prereqs (NOT redistributed — see data/reference/README.md):
  - MAGMA binary: on PATH, or data/reference/magma(.exe), or env MAGMA=...
  - data/reference/NCBI37.3.gene.loc        (gene locations, GRCh37 — matches GWAS)
  - data/reference/g1000_eur.{bed,bim,fam}   (1000G EUR reference panel, hg19)

GWAS build is GRCh37/hg19 (verified) → no liftOver needed.

Usage:
  python scripts/01_gene_based_magma.py                 # primary pheno from config
  python scripts/01_gene_based_magma.py --pheno alphaOcc
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd

from utils import ROOT, load_config, log_run

MAGMA_DIR = ROOT / "results" / "tables" / "magma"
REFERENCE = ROOT / "data" / "reference"


def find_magma() -> str:
    """Locate the MAGMA binary: env MAGMA, then data/reference/, then PATH."""
    env = os.environ.get("MAGMA")
    if env and Path(env).exists():
        return env
    for name in ("magma.exe", "magma"):
        cand = REFERENCE / name
        if cand.exists():
            return str(cand)
    onpath = shutil.which("magma")
    if onpath:
        return onpath
    raise SystemExit(
        "MAGMA binary not found. Download the Windows build from "
        "https://ctg.cncr.nl/software/magma and place magma.exe in data/reference/ "
        "(or set env MAGMA=...). See data/reference/README.md."
    )


def check_reference() -> tuple[Path, Path]:
    gene_loc = REFERENCE / "NCBI37.3.gene.loc"
    panel = REFERENCE / "g1000_eur"
    missing = [str(p) for p in (gene_loc,) if not p.exists()]
    if not (REFERENCE / "g1000_eur.bed").exists():
        missing.append(str(panel) + ".{bed,bim,fam}")
    if missing:
        raise SystemExit(
            "Missing MAGMA reference files:\n  " + "\n  ".join(missing) +
            "\nDownload from https://ctg.cncr.nl/software/magma "
            "(gene locations NCBI37.3, and 1000G European panel). "
            "See data/reference/README.md."
        )
    return gene_loc, panel


def write_snp_loc(sumstats: Path, cols: dict, out: Path) -> int:
    """Write a MAGMA snp-loc file: SNP CHR BP (whitespace, no header)."""
    usecols = [cols["snp"], cols["chr"], cols["bp"]]
    df = pd.read_csv(sumstats, sep=r"\s+", usecols=usecols)
    df = df[usecols]  # enforce SNP, CHR, BP order
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, sep="\t", header=False, index=False)
    return len(df)


def run(cmd: list[str]) -> None:
    print("[run]", " ".join(str(c) for c in cmd))
    subprocess.run(cmd, check=True)


def main() -> None:
    cfg = load_config()
    g = cfg["gwas"]
    ap = argparse.ArgumentParser()
    ap.add_argument("--pheno", default=g["primary"], choices=list(g["files"]))
    args = ap.parse_args()

    sumstats = ROOT / g["files"][args.pheno]
    if not sumstats.exists():
        raise SystemExit(f"Sumstats not found: {sumstats}")

    magma = find_magma()
    gene_loc, panel = check_reference()
    cols = g["cols"]
    window = ",".join(str(w) for w in cfg["magma"]["window_kb"])
    prefix = MAGMA_DIR / args.pheno

    # 1. SNP location file
    snploc = MAGMA_DIR / f"{args.pheno}.snploc"
    n_snps = write_snp_loc(sumstats, cols, snploc)
    print(f"[snploc] {snploc}  ({n_snps} SNPs)")

    # 2. Annotate SNPs -> genes
    run([magma, "--annotate", f"window={window}",
         "--snp-loc", str(snploc), "--gene-loc", str(gene_loc),
         "--out", str(prefix)])

    # 3. Gene-based test (SNP-wise, per-SNP N)
    run([magma, "--bfile", str(panel),
         "--gene-annot", f"{prefix}.genes.annot",
         "--pval", str(sumstats), f"use={cols['snp']},{cols['pval']}", f"ncol={cols['n']}",
         "--out", str(prefix)])

    print(f"[done] gene-based results: {prefix}.genes.out")
    print("Next: python scripts/01_postprocess_magma.py --pheno", args.pheno)

    log_run("01_gene_based_magma", {
        "phenotype": args.pheno, "sumstats": g["files"][args.pheno],
        "n_snps": n_snps, "window_kb": cfg["magma"]["window_kb"],
        "build": g["build"], "magma": magma})


if __name__ == "__main__":
    main()
