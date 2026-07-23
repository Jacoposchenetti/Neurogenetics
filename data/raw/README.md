# Raw external data — NOT redistributed (download yourself)

These datasets have their own licenses / access agreements. This folder is
git-ignored except for this file.

## 1. ENIGMA-EEG GWAS summary statistics (Smit et al. 2018) — Fase 1  [OBTAINED]

- Requested via the Google Form linked at
  https://enigma.ini.usc.edu/research/download-enigma-gwas-results/ (reviewed ~5 business days),
  then downloaded from the approved `website_downloads/ENIGMA-EEG/` directory.
- 6 phenotypes, 20181101 release, `select_N6000_Hom100` (per-SNP N ≈ 6000–7700):
  alphaCz, alphaOcc, peakOcc (alpha peak freq), thetaCz, betaCz, deltaCz.
- Files live one-per-folder under `data/raw/ENIGMA-EEG_20181101_<pheno>.txt/select_N6000_Hom100.<pheno>.txt`.
  Paths are recorded in `config.yaml:gwas.files`.
- **Columns (13):** `RS A1 A2 N Z PVAL DIR HomIsq HomSE HomDf HomP CHR BP`.
  MAGMA mapping: SNP=RS, CHR, BP, P=PVAL, N (per-SNP).
- **Genome build = GRCh37/hg19 (VERIFIED:** rs3094315 @ chr1:752566). No liftOver —
  matches the NCBI37 gene.loc + g1000_eur panel.
- Restricted-access data: NOT redistributed (folder is git-ignored).

## 2. Allen Human Brain Atlas — Fase 3

- **Do not download manually.** `abagen.fetch_microarray()` / `get_expression_data()`
  downloads and caches it automatically on first run (~4 GB, 6 donors).
- Cache location is git-ignored (`abagen-data/`).

## 3. (Optional) HCP-MEG — secondary heritability analysis (Fase 5.3)

- Requires the HCP Data Use Agreement (free) at https://db.humanconnectome.org.
