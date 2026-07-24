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

Run Phase 1 with: `python scripts/01_gene_based_magma.py --pheno alphaCz`.
Exact versions/params are logged to `results/logs/` on each run.

Note: only `scripts/01_*` need this folder. The rest of the pipeline (expression,
enrichment, figures) runs from the committed derived results.

## Yan-600 homotopic parcellation (Phase 6)

Public download, git-ignored.

Needed only by `scripts/06_parcellation_sensitivity.py --parc yan600`
(`schaefer100` is fetched automatically by `netneurotools`; no action required).

Source: Yan et al. 2023 homotopic parcellation in the CBIG repository, the
fsaverage5 / Kong2022-17-networks variant. **The scripts expect different
filenames than the source**, so rename on download:

| Download from CBIG | Save as |
|---|---|
| `lh.600Parcels_Kong2022_17Networks.annot` | `data/reference/yan600_lh.annot` |
| `rh.600Parcels_Kong2022_17Networks.annot` | `data/reference/yan600_rh.annot` |

```bash
CBIG=https://raw.githubusercontent.com/ThomasYeoLab/CBIG/master/stable_projects/brain_parcellation/Yan2023_homotopic/parcellations/FreeSurfer/fsaverage5/label/kong17
curl -L $CBIG/lh.600Parcels_Kong2022_17Networks.annot -o data/reference/yan600_lh.annot
curl -L $CBIG/rh.600Parcels_Kong2022_17Networks.annot -o data/reference/yan600_rh.annot
```

Sanity check: each file is an fsaverage5 annot with 10,242 vertices and 300
parcels per hemisphere (600 total).
