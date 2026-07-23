# Neurogenetics: imaging transcriptomics of the cortical alpha-rhythm sources targeted by closed-loop EEG–TMS

Code and derived results for the study:

> **A significant enrichment that is not: spatial nulls, co-expression, and the
> imaging transcriptomics of EEG oscillatory-power genetics**

The manuscript is in [`paper/manuscript.pdf`](paper/manuscript.pdf) (source:
`paper/manuscript.md`, `paper/manuscript.tex`).

## What this study does

It asks two questions at once.

**Substantive:** are genes associated with EEG alpha power expressed more highly in
the 41 cortical generators of the alpha rhythm (Tabarelli et al., 2022) than in the
rest of cortex, and is any such enrichment *specific* to the alpha band?

**Methodological:** does the answer survive complementary null models — or only the
spin test that is the field standard?

Pipeline: ENIGMA-EEG GWAS summary statistics → MAGMA gene-based analysis → Allen
Human Brain Atlas regional expression (`abagen`) → enrichment tested against a
spatial (spin) null, a co-expression-aware gene-set null, and a positive control,
across 6 phenotypes and 3 parcellations.

**Result.** Judged by the spin test alone the hypothesis appears supported
(*p*_spin = 0.022), and the result is corroborated by an independent surrogate model
while the pipeline passes a positive control. It then dissolves: the enrichment is
not band-specific, does not replicate across alpha phenotypes, does not survive a
change of parcellation or FDR correction, and — decisively — is indistinguishable
from random gene sets of matched size (*p*_geneset = 0.33). Correcting for spatial
autocorrelation is necessary but not sufficient; enrichment claims should report a
gene-set null alongside the spatial null.

---

## Repository layout

```
config.yaml              single source of all analysis parameters (seeds, atlases, thresholds)
requirements.txt         pinned Python environment
scripts/                 analysis pipeline (numbered in execution order)
data/raw/                GWAS summary statistics  (NOT redistributed - see below)
data/reference/          MAGMA binary, 1000G panel, atlases (NOT redistributed)
data/derived/            gene-level statistics and gene sets produced by the pipeline
results/tables/          MAGMA output, enrichment statistics
results/figures/         intermediate figures produced by the analysis scripts
results/logs/            per-run JSON logs (parameters, seeds, package versions)
paper/                   manuscript + publication figures
paper/figures/           final figures and shared plotting style
docs/methods_notes.md    methodological decisions and rationale
```

## Requirements

- **Python 3.12.** The neuroimaging stack does not provide wheels for 3.13+; 3.12 is
  required. Verified on Windows 11, but nothing is platform-specific except the
  MAGMA binary.
- **MAGMA v1.10** (external C binary, not pip-installable).
- ~10 GB free disk (Allen Human Brain Atlas cache ~6.4 GB, 1000G reference panel ~3.6 GB).

```bash
python3.12 -m venv venv
# Windows:  venv\Scripts\activate
# Linux/macOS:  source venv/bin/activate
pip install -r requirements.txt
```

Two notes on the environment, both already reflected in `requirements.txt`:

- `abagen` is installed from GitHub, not PyPI: the released 0.1.3 calls
  `DataFrame.set_axis(..., inplace=)`, removed in pandas 2, and crashes.
- `setuptools<81` is pinned because `abagen` still imports `pkg_resources`.

---

## Data you must obtain yourself

Raw data are **not redistributed** here — they carry their own licences and access
agreements. Only derived results are versioned. Each folder has a README with the
exact instructions.

| Data | Source | Access | Goes in |
|---|---|---|---|
| ENIGMA-EEG GWAS summary statistics | ENIGMA-EEG working group | application form, reviewed ~5 business days | `data/raw/` |
| MAGMA v1.10 binary + NCBI37.3 gene locations + 1000G EUR panel | CTG lab (`cncr.nl/research/magma/`) | public download | `data/reference/` |
| Allen Human Brain Atlas microarray | Allen Institute | downloaded automatically by `abagen` on first run (~6.4 GB, cached) | (cache) |
| Glasser HCP-MMP1.0, Schaefer-100 | fetched automatically by `netneurotools` | public | (cache) |
| Yan-600 homotopic parcellation | CBIG repository (`ThomasYeoLab/CBIG`) | public download | `data/reference/` |

Expected filenames are documented in [`data/raw/README.md`](data/raw/README.md) and
[`data/reference/README.md`](data/reference/README.md).

**Genome build:** the ENIGMA-EEG statistics are GRCh37/hg19 (verified: rs3094315 at
chr1:752,566), matching the Build-37 gene locations and the 1000G panel — no
liftOver is required. If you substitute another GWAS, check this first.

---

## Reproducing the analysis

All parameters live in `config.yaml`; the scripts read from it, so the pipeline is
reproduced by running them in order. Every script writes a JSON log to
`results/logs/` recording its parameters, the random seed and installed package
versions.

```bash
# 1. Gene-based association (per phenotype). Requires MAGMA + reference panel.
#    ~80 min per phenotype; run for all six.
python scripts/01_gene_based_magma.py   --pheno alphaCz
python scripts/01_postprocess_magma.py  --pheno alphaCz     # -> gene sets, gene-level Z
#    repeat for: alphaOcc  peakOcc  thetaCz  betaCz  deltaCz

# 2. Target region set: the 41 Tabarelli alpha-source ROIs in Glasser space
python scripts/02_roi_mapping.py

# 2b. Resample the Glasser atlas from fsLR32k to fsaverage5 (abagen's sample space)
python scripts/make_glasser_fsaverage5.py

# 3. Regional gene expression from the Allen Human Brain Atlas
#    First run downloads the AHBA (~6.4 GB, cached afterwards).
python scripts/03_abagen_expression.py --atlas glasser

# 4. Primary spatial enrichment test (spin null, 10,000 rotations)
python scripts/04_spin_test_enrichment.py --pheno alphaCz

# 5. Specificity across all phenotypes + FDR over the full grid
python scripts/05_sensitivity.py

# 6. Parcellation sensitivity
python scripts/06_parcellation_sensitivity.py --parc schaefer100 --pheno alphaCz
python scripts/06_parcellation_sensitivity.py --parc yan600      --pheno alphaCz

# 7. Method-validation analyses (positive control, gene-set null, surrogate check)
python scripts/07_method_validation.py

# 8. Figures
python scripts/08_figures_brain.py     # surface figures  (Fig 2, Fig 5)
python scripts/09_figures_stats.py     # statistical figures (Fig 1, 3, 4, 6)
```

`scripts/enrichment_lib.py` holds the shared machinery (parcel centroids, spin
generation, the two alpha scores, the spin p-value) used by steps 4–7; `scripts/utils.py`
handles config loading and run logging.

### Determinism

A single seed (`config.yaml: seed`) governs every permutation. Spin permutations are
cached to `results/derived/` and reused, so repeated runs give identical numbers.
Re-running the pipeline from scratch reproduces every value in the manuscript
tables and figures.

---

## Rebuilding the manuscript

The manuscript is written in Markdown and converted to LaTeX, then compiled to PDF:

```bash
cd paper
pandoc manuscript.md -s -V geometry:margin=1in -V fontsize=11pt -o manuscript.tex
tectonic -X compile manuscript.tex
```

Requires [pandoc](https://pandoc.org) and [tectonic](https://tectonic-typesetting.github.io).
Figure paths in `manuscript.md` are relative to `paper/`, so run the commands from
inside that directory.

---

## Key methodological points

- **Parcellation.** Glasser HCP-MMP1.0 is the primary atlas because the alpha sources
  are natively defined in it, making the region mapping one-to-one and avoiding a
  subjective region translation. Schaefer-100 and Yan-600 are sensitivity analyses,
  with the target set re-mapped by fsaverage5 surface overlap (≥50% of a parcel's
  cortical vertices).
- **Two enrichment scores.** A threshold-free *continuous* score (rank correlation
  between a region's expression profile and the genome-wide gene-level statistics)
  is the primary readout, since only one gene survives genome-wide correction; a
  *top-100* gene-set composite is the sensitivity analysis.
- **Spatial null.** Enrichment is tested against sphere-rotation (spin) permutations
  of parcel centroids on fsaverage5, implemented in `enrichment_lib.gen_spins` because
  the spin generator was removed from the installed `netneurotools` version. It is
  cross-checked against `brainsmash` surrogate maps.
- **Co-expression null.** Because transcriptomic enrichment is prone to false
  positives from within-set co-expression, the top-100 result is additionally tested
  against 10,000 random gene sets of matched size.
- **AHBA donors.** Five of six donors are used; donor 15496 was unavailable at source
  (the Allen file returns HTTP 404). Set `config.yaml: abagen.donors` back to all six
  if it is restored.

Further rationale is in [`docs/methods_notes.md`](docs/methods_notes.md).

## Citation

If you use this code, please cite the manuscript and the underlying resources:
ENIGMA-EEG (Smit et al., 2018), the Allen Human Brain Atlas (Hawrylycz et al., 2012),
`abagen` (Markello et al., 2021), MAGMA (de Leeuw et al., 2015), and the parcellations
(Glasser et al., 2016; Schaefer et al., 2018; Yan et al., 2023). Full references are
in the manuscript.

## Licence

Code is released under the MIT Licence (see `LICENSE`). Raw GWAS, AHBA and reference
panel data remain under their original terms and are not redistributed here.
