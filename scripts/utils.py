"""Shared helpers: config loading + reproducibility logging.

Every pipeline script imports these so parameters and package versions are
captured in results/logs/ for full reproducibility.
"""
from __future__ import annotations

import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = ROOT / "config.yaml"
LOG_DIR = ROOT / "results" / "logs"


def load_config(path: Path | str = CONFIG_PATH) -> dict:
    with open(path, "r", encoding="utf-8") as fh:
        return yaml.safe_load(fh)


def _pkg_versions() -> dict:
    """Best-effort version capture of the scientific stack."""
    versions = {"python": sys.version.split()[0], "platform": platform.platform()}
    for name in ("numpy", "pandas", "scipy", "statsmodels", "nibabel",
                 "nilearn", "abagen", "netneurotools", "brainsmash"):
        try:
            mod = __import__(name)
            versions[name] = getattr(mod, "__version__", "unknown")
        except Exception:
            versions[name] = "not installed"
    return versions


def log_run(phase: str, params: dict) -> Path:
    """Write a timestamped JSON log of parameters + environment for one run."""
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    record = {
        "phase": phase,
        "timestamp_utc": stamp,
        "params": params,
        "environment": _pkg_versions(),
    }
    out = LOG_DIR / f"{phase}_{stamp}.json"
    with open(out, "w", encoding="utf-8") as fh:
        json.dump(record, fh, indent=2)
    print(f"[log] {out}")
    return out
