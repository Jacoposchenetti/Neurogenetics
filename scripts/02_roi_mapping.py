"""Fase 2 — Map Tabarelli et al. (2022) alpha-source ROIs onto the parcellation.

The paper is natively defined in the Glasser HCP-MMP1.0 atlas (360 regions), so
the "region-in-paper -> standard parcel" mapping is essentially 1:1 by label.
This script produces a binary target / non-target vector over the parcellation
and writes it to data/derived/ for use by Fase 3-4.

Run:  python scripts/02_roi_mapping.py
Output: data/derived/roi_target_glasser.csv  (columns: label, sector, is_target, is_hub)
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from utils import ROOT, load_config, log_run

# --- 41 alpha-source ROIs, HCP-MMP1.0 labels, from Tabarelli et al. 2022 Table 1 ---
# Grouped by the paper's three cortical sectors. Bilateral (applied to L and R).
ROI_BY_SECTOR: dict[str, list[str]] = {
    "occipital": [
        "V1", "V2", "ProS", "V3", "V4", "V6", "V6A", "V7",
        "IPS1", "V3A", "V3B", "V3CD", "IP0", "PGp", "LO1", "LO2",
    ],
    "parietal": ["1", "2", "3a", "3b", "4", "6mp", "6d"],
    "frontal": [
        "8BL", "9p", "9m", "9a", "8Ad", "9-46d", "8BM", "8Av", "46", "8C",
        "p9-46v", "a32pr", "d32", "a9-46v", "10d", "p10p", "p47r", "IFSa",
    ],
}

# Prefrontal connectivity hubs (targeted H1 test). See config.yaml:roi.hub_subset.
HUB_LABELS = {"8C", "8Av", "9a"}


def build_roi_table() -> pd.DataFrame:
    rows = []
    for sector, labels in ROI_BY_SECTOR.items():
        for lab in labels:
            rows.append(
                {
                    "label": lab,
                    "sector": sector,
                    "is_target": 1,
                    "is_hub": int(lab in HUB_LABELS),
                }
            )
    df = pd.DataFrame(rows)
    assert len(df) == 41, f"expected 41 ROIs, got {len(df)}"
    return df


# NOTE (Fase 3 hand-off): abagen returns expression indexed by parcel id, not by
# Glasser label. The label<->parcel-id crosswalk depends on the atlas image you
# feed abagen (e.g. the HCP-MMP volumetric NIfTI or fsaverage annot). Resolve it
# THERE, using the atlas' own label table, then left-join this table on `label`.
# Glasser labels here are hemisphere-agnostic (e.g. "V1"): the atlas typically
# encodes L_V1 / R_V1 — expand bilaterally when joining unless a lateralized
# hypothesis is intended.


def main() -> None:
    cfg = load_config()
    df = build_roi_table()
    out = ROOT / "data" / "derived" / "roi_target_glasser.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)
    print(f"[write] {out}  ({len(df)} ROIs; {df.is_hub.sum()} hubs)")
    print(df.groupby("sector").size().to_string())

    log_run(
        "02_roi_mapping",
        {
            "n_roi": int(len(df)),
            "n_hub": int(df.is_hub.sum()),
            "parcellation": cfg["parcellation"]["primary"],
            "source": "Tabarelli et al. 2022, Brain Sci 12(3):348, Table 1",
            "output": str(out.relative_to(ROOT)),
        },
    )


if __name__ == "__main__":
    main()
