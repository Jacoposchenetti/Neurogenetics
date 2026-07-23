# Reference panels & external binaries — NOT redistributed

Git-ignored except for this file.

## MAGMA tooling — Fase 1  [OBTAINED, git-ignored]

Downloaded 2026-07-21 from the CTG lab (https://cncr.nl/research/magma/),
MAGMA **v1.10**. GWAS build is GRCh37/hg19, so the Build-37 gene locations +
European panel match (no liftOver).

Present in this folder (all git-ignored):
```
magma.exe                 MAGMA v1.10 Windows 64-bit static
magma_manual_v1.10.pdf    reference manual
NCBI37.3.gene.loc         gene locations, build 37
g1000_eur.bed/.bim/.fam   1000 Genomes European panel (503 indiv, 22.6M SNPs)
g1000_eur.synonyms        SNP synonyms (auto-detected by MAGMA)
```

Source URLs (MAGMA v1.10, SURFdrive shares; append `/download`):
- Windows binary: https://vu.data.surfsara.nl/index.php/s/TOH4SuvczAKE29d
- Gene loc Build 37: https://vu.data.surfsara.nl/index.php/s/Pj2orwuF2JYyKxq
- European reference: https://vu.data.surfsara.nl/index.php/s/VZNByNwpD8qqINe

Run Fase 1 with: `python scripts/01_gene_based_magma.py --pheno alphaCz`.
Exact versions/params are logged to `results/logs/` on each run.
