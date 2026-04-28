"""
utils.py — Shared helpers: file validation, threshold logic, config loading.
"""

from __future__ import annotations

import json
import os
import shutil
import sys
from pathlib import Path
from typing import Any


# ---------------------------------------------------------------------------
# File / tool validation
# ---------------------------------------------------------------------------

def check_file_readable(path: str, label: str = "File") -> Path:
    """Assert that *path* exists and is readable; return a Path object."""
    p = Path(path)
    if not p.exists():
        _fatal(f"{label} not found: {path}")
    if not p.is_file():
        _fatal(f"{label} is not a regular file: {path}")
    if not os.access(p, os.R_OK):
        _fatal(f"{label} is not readable: {path}")
    return p


def check_bam_index(bam_path: str) -> None:
    """Ensure a BAM index (.bai) exists next to the BAM file.

    Checks two common index naming conventions:
      sample.bam  →  sample.bam.bai   (samtools default)
      sample.bam  →  sample.bai       (alternative)
    """
    bam = Path(bam_path)
    # .with_suffix replaces only the last extension, so for "sample.bam"
    # we get "sample.bai" — that's the alternative form. Append ".bai"
    # directly to get "sample.bam.bai" (samtools default).
    bai_same = Path(str(bam) + ".bai")   # sample.bam.bai  (samtools default)
    bai_alt  = bam.with_suffix(".bai")   # sample.bai       (alternative)
    if not bai_same.exists() and not bai_alt.exists():
        _fatal(
            f"BAM index not found for {bam_path}.\n"
            f"Run:  samtools index {bam_path}"
        )


def check_mosdepth(sif: str | None = None) -> tuple[str, str | None]:
    """
    Return (mosdepth_exe_or_cmd, sif_path_or_None).

    Resolution order:
      1. If --mosdepth-sif PATH given → use Singularity
      2. If mosdepth in PATH → use directly
      3. Auto-detect well-known SIF locations
      4. Fatal error with install hint
    """
    # 1. Explicit SIF path
    if sif:
        sif_path = Path(sif)
        if not sif_path.exists():
            _fatal(f"mosdepth SIF image not found: {sif}")
        singularity = shutil.which("singularity") or shutil.which("apptainer")
        if singularity is None:
            _fatal(
                "singularity/apptainer not found in PATH but --mosdepth-sif was given.\n"
                "Install Singularity or Apptainer, or put mosdepth directly in PATH."
            )
        return singularity, str(sif_path)

    # 2. mosdepth in PATH
    exe = shutil.which("mosdepth")
    if exe is not None:
        return exe, None

    # 3. Auto-detect well-known SIF locations (generic HPC / workstation paths)
    _WELL_KNOWN_SIFS = [
        os.path.expanduser("~/singularity/mosdepth.sif"),
        "/opt/software/mosdepth/mosdepth.sif",
        "/usr/local/lib/mosdepth/mosdepth.sif",
        "/opt/mosdepth/mosdepth.sif",
    ]
    for candidate in _WELL_KNOWN_SIFS:
        if Path(candidate).exists():
            singularity = shutil.which("singularity") or shutil.which("apptainer")
            if singularity:
                print(
                    f"[vcfdash] mosdepth not in PATH — using SIF: {candidate}",
                    file=__import__("sys").stderr,
                )
                return singularity, candidate

    # 4. Fatal
    _fatal(
        "mosdepth not found in PATH and no SIF image located.\n"
        "Options:\n"
        "  a) conda install -c bioconda mosdepth\n"
        "  b) vcfdash --mosdepth-sif /path/to/mosdepth.sif ..."
    )
    return "", None  # unreachable


def check_output_dir(outdir: str) -> Path:
    """Create output directory if it does not exist; return Path."""
    p = Path(outdir)
    p.mkdir(parents=True, exist_ok=True)
    return p


# ---------------------------------------------------------------------------
# Config file
# ---------------------------------------------------------------------------

DEFAULT_CONFIG: dict[str, Any] = {
    "min_dp":  20,
    "min_dp2": 10,
    "pass_pct_20x": 95.0,
    "warn_pct_20x": 80.0,
    "pass_depth_fraction": 1.0,
    "warn_depth_fraction": 0.5,
}


def load_config(config_path: str | None, args) -> dict[str, Any]:
    """
    Merge CLI args with optional JSON config file.
    CLI args always take precedence over config file values.
    """
    cfg = dict(DEFAULT_CONFIG)

    if config_path:
        config_path_p = check_file_readable(config_path, "Config file")
        try:
            with open(config_path_p) as fh:
                file_cfg = json.load(fh)
            cfg.update(file_cfg)
        except json.JSONDecodeError as exc:
            _fatal(f"Config file is not valid JSON: {exc}")

    # CLI args always win
    cfg["min_dp"]  = args.min_dp
    cfg["min_dp2"] = args.min_dp2

    return cfg


# ---------------------------------------------------------------------------
# Threshold / status logic
# ---------------------------------------------------------------------------

def assign_gene_status(
    mean_depth: float,
    pct_20x: float,
    min_dp: int,
    cfg: dict[str, Any] | None = None,
) -> str:
    """
    Return 'PASS', 'WARN', or 'FAIL' for a gene coverage region.

    PASS : mean_depth >= min_dp  AND pct_20x >= 95
    WARN : mean_depth >= min_dp * 0.5  AND pct_20x >= 80
    FAIL : anything below WARN thresholds
    """
    if cfg is None:
        cfg = DEFAULT_CONFIG

    pass_pct   = cfg.get("pass_pct_20x", 95.0)
    warn_pct   = cfg.get("warn_pct_20x", 80.0)
    pass_frac  = cfg.get("pass_depth_fraction", 1.0)
    warn_frac  = cfg.get("warn_depth_fraction", 0.5)

    if mean_depth >= min_dp * pass_frac and pct_20x >= pass_pct:
        return "PASS"
    if mean_depth >= min_dp * warn_frac and pct_20x >= warn_pct:
        return "WARN"
    return "FAIL"


def coverage_color_class(status: str) -> str:
    """Map status string to a CSS class name."""
    return {"PASS": "status-pass", "WARN": "status-warn", "FAIL": "status-fail"}.get(
        status, "status-fail"
    )


def metric_color_class(value: float, good: float, warn: float, higher_is_better: bool = True) -> str:
    """Return a CSS class based on whether a numeric metric is good/warn/bad."""
    if higher_is_better:
        if value >= good:
            return "metric-good"
        if value >= warn:
            return "metric-warn"
        return "metric-bad"
    else:
        if value <= good:
            return "metric-good"
        if value <= warn:
            return "metric-warn"
        return "metric-bad"


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _fatal(message: str) -> None:
    """Print an error message to stderr and exit(1)."""
    print(f"ERROR: {message}", file=sys.stderr)
    sys.exit(1)


def safe_float(value: Any, default: float = 0.0) -> float:
    """Safely convert *value* to float, returning *default* on failure."""
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def safe_int(value: Any, default: int = 0) -> int:
    """Safely convert *value* to int, returning *default* on failure."""
    try:
        return int(value)
    except (TypeError, ValueError):
        return default
