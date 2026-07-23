# Methods notes (for the paper's Methods + reproducibility)

## Parcellation & ROI definition (Fase 2)

- Primary atlas: **Glasser HCP-MMP1.0** (360 cortical parcels, 180/hemisphere),
  matching Tabarelli et al. 2022, which defines the alpha sources natively in this
  space. This makes the "paper-region → parcel" mapping ~1:1 by label and removes
  the subjective-mapping concern.
- **41 alpha-source ROIs** taken from Tabarelli et al. 2022, Table 1
  (16 occipital, 7 parietal, 18 frontal). All 41 labels are present in the atlas;
  each maps bilaterally (L+R) → 82 parcels flagged as target.
- Prefrontal connectivity hubs (8C, 8Av, 9a) additionally flagged for a targeted
  H1 test, since the paper's connectivity is hub-weighted, not uniform.

## Glasser fsLR32k → fsaverage5 resampling (make_glasser_fsaverage5.py)

abagen maps AHBA microarray samples to **fsaverage5**, but Glasser is native to
**fsLR32k**. The two use different spherical registrations, so a naive sphere
nearest-neighbor would be spatially wrong.

We resample using the canonical HCP `resample_fsaverage` sphere pair (from
`netneurotools.fetch_hcp_standards`, i.e. the WU standard_mesh_atlases), in which
both surfaces are expressed in the **common fsaverage spherical frame**:

- `fs_LR-deformed_to-fsaverage.{L,R}.sphere.32k_fs_LR` — the 32 492 fsLR vertices
  (same order as the mmpall Glasser labels) warped into the fsaverage frame.
- `fsaverage5_std_sphere.{L,R}.10k_fsavg` — the 10 242 fsaverage5 vertices (same
  order as FreeSurfer fsaverage5, which abagen expects) in the same frame.

For each fsaverage5 vertex we assign the label of the nearest fsLR vertex
(scipy cKDTree, k=1) in that shared frame. This mirrors what Connectome Workbench
does, without requiring the native `wb_command` binary. Nearest-neighbor is the
standard choice for resampling discrete label maps.

**Validation:** both hemispheres yield exactly **180 parcels** post-resampling
(no parcel lost), with ~9 370 cortical vertices/hemi and ~870 medial-wall vertices
— consistent with fsaverage5. Trusted sources only (OSF via netneurotools).

## AHBA expression (Fase 3, abagen)

- abagen params logged in `config.yaml:abagen` and every run's `results/logs/`.
- **Donor exclusion:** the analysis uses **5 of 6** AHBA donors
  (9861, 10021, 12876, 14380, 15697). Donor **15496 (H0351.1015)** is excluded
  because its Allen normalized-microarray file (well_known_file 178238266) returns
  **HTTP 404** at the source as of 2026-07-21 — the Allen DB record exists but the
  download route is broken server-side; no alternative URL is published. 15496 is
  a left-hemisphere-only donor. **Action:** re-run with all 6 donors once Allen
  restores the file (flip `config.yaml:abagen.donors` back to all six).

## Reproducibility

- Global seed, parcellation, abagen params, spin-test settings: `config.yaml`.
- Every script writes a JSON run-log (params + package versions) to `results/logs/`.
- Raw data (AHBA, ENIGMA-EEG, 1000G) are not redistributed; download instructions
  live in `data/raw/README.md` and `data/reference/README.md`.
