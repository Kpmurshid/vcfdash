"""
coverage.py — mosdepth wrapper, per-base and per-region coverage parsing.
"""

from __future__ import annotations

import gzip
import os
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any

from .utils import assign_gene_status, safe_float, safe_int, _fatal, check_mosdepth


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def run_mosdepth(
    bam: str,
    bed: str,
    outdir: str,
    prefix: str,
    min_dp: int = 20,
    min_dp2: int = 10,
    no_sparklines: bool = False,
    mosdepth_sif: str | None = None,
) -> dict[str, Any]:
    """
    Run mosdepth on *bam* restricted to *bed* regions.

    If *mosdepth_sif* is given (or auto-detected), mosdepth is run via
    ``singularity exec <sif> mosdepth …``.  Otherwise mosdepth must be in PATH.

    Returns a dict with keys:
      - per_region   : list of region coverage records
      - global_stats : dict with mean_coverage, pct_above_20x, pct_above_10x
      - sparkline_data : dict  chrom:pos -> list of (position, depth) tuples
                         (empty when no_sparklines=True)
    """
    out_prefix = str(Path(outdir) / prefix)
    quantize_arg = f"0:{min_dp2}:{min_dp}:"

    # Resolve mosdepth executable / SIF
    singularity_exe, resolved_sif = check_mosdepth(mosdepth_sif)

    # Build the mosdepth argument list
    # mosdepth 0.3.x flags:
    #   -b/--by          BED file for per-region output
    #   -q/--quantize    quantize output (segments separated by colons)
    #   -n/--no-per-base skip per-base depth (fast mode, no sparklines)
    #   -x/--fast-mode   skip cigar ops / mate overlap correction (faster)
    mosdepth_args = [
        "--by", bed,
        "--quantize", quantize_arg,
        "--fast-mode",          # always use fast mode for speed
    ]
    if no_sparklines:
        mosdepth_args += ["--no-per-base"]   # skip per-base.bed.gz entirely
    mosdepth_args += [out_prefix, bam]

    # Prefix with singularity exec if using a SIF image
    if resolved_sif:
        # Bind the directories containing input files so Singularity can read them
        bam_dir    = str(Path(bam).resolve().parent)
        bed_dir    = str(Path(bed).resolve().parent)
        outdir_abs = str(Path(outdir).resolve())
        binds      = ",".join({bam_dir, bed_dir, outdir_abs})
        # Probe the SIF for the mosdepth binary location (may be /opt/mosdepth or
        # in PATH as 'mosdepth'). Try /opt/mosdepth first, then fall back to 'mosdepth'.
        mosdepth_in_sif = _find_mosdepth_in_sif(singularity_exe, resolved_sif)
        cmd = [
            singularity_exe, "exec",
            "--bind", binds,
            resolved_sif,
            mosdepth_in_sif,
        ] + mosdepth_args
    else:
        cmd = [singularity_exe] + mosdepth_args  # singularity_exe IS mosdepth here

    print(f"[vcfdash] Running mosdepth …", file=sys.stderr)
    if resolved_sif:
        print(f"[vcfdash]   via Singularity SIF: {resolved_sif}", file=sys.stderr)

    try:
        result = subprocess.run(
            cmd,
            check=True,
            capture_output=True,
            text=True,
        )
        if result.stderr:
            print(result.stderr, file=sys.stderr)
    except subprocess.CalledProcessError as exc:
        _fatal(
            f"mosdepth failed (exit code {exc.returncode}).\n"
            f"stderr:\n{exc.stderr}"
        )

    print(f"[vcfdash] Parsing mosdepth output …", file=sys.stderr)
    return parse_mosdepth_output(
        prefix=out_prefix,
        min_dp=min_dp,
        min_dp2=min_dp2,
        no_sparklines=no_sparklines,
    )


def parse_mosdepth_output(
    prefix: str,
    min_dp: int = 20,
    min_dp2: int = 10,
    no_sparklines: bool = False,
) -> dict[str, Any]:
    """
    Parse mosdepth output files located at *prefix*.*.

    Expected files:
      {prefix}.regions.bed.gz        — per-region mean depths
      {prefix}.mosdepth.summary.txt  — global stats
      {prefix}.per-base.bed.gz       — per-base depth (optional, for sparklines)
    """
    regions_file  = Path(f"{prefix}.regions.bed.gz")
    summary_file  = Path(f"{prefix}.mosdepth.summary.txt")
    perbase_file  = Path(f"{prefix}.per-base.bed.gz")

    per_region   = _parse_regions_bed(regions_file, min_dp, min_dp2)
    # If regions have no real gene names (coordinate-only labels), aggregate
    # them by chromosome to keep the report small and usable.
    per_region   = _aggregate_regions(per_region, min_dp, min_dp2)
    global_stats = _parse_summary(summary_file, min_dp, min_dp2, per_region)

    sparkline_data: dict[str, list] = {}
    if not no_sparklines and perbase_file.exists():
        sparkline_data = _build_sparkline_index(perbase_file)

    return {
        "per_region":    per_region,
        "global_stats":  global_stats,
        "sparkline_data": sparkline_data,
    }


# ---------------------------------------------------------------------------
# Region BED parsing
# ---------------------------------------------------------------------------

def _parse_regions_bed(
    path: Path,
    min_dp: int,
    min_dp2: int,
) -> list[dict[str, Any]]:
    """Parse {prefix}.regions.bed.gz → list of per-region dicts."""
    if not path.exists():
        _fatal(f"mosdepth regions file not found: {path}")

    records: list[dict[str, Any]] = []

    with gzip.open(path, "rt") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue

            chrom = parts[0]
            start = safe_int(parts[1])
            end   = safe_int(parts[2])

            # Column 3 may be gene name (4-col BED) or mean depth (mosdepth output)
            # mosdepth outputs: chrom start end name mean_depth  (5 cols) when BED has name
            # or:               chrom start end mean_depth       (4 cols) when no name
            if len(parts) >= 5:
                gene       = parts[3] if parts[3] else f"{chrom}:{start}-{end}"
                mean_depth = safe_float(parts[4])
            else:
                gene       = f"{chrom}:{start}-{end}"
                mean_depth = safe_float(parts[3])

            # pct_20x and pct_10x are not directly in the regions bed; we compute
            # them later from per-base data if available, or estimate from summary.
            # For now, store placeholders — they will be updated by _enrich_regions().
            records.append({
                "chrom":      chrom,
                "start":      start,
                "end":        end,
                "gene":       gene,
                "mean_depth": round(mean_depth, 2),
                "pct_20x":    None,   # populated later
                "pct_10x":    None,
                "status":     None,   # populated later
            })

    return records


# ---------------------------------------------------------------------------
# Summary file parsing
# ---------------------------------------------------------------------------

def _parse_summary(
    path: Path,
    min_dp: int,
    min_dp2: int,
    per_region: list[dict],
) -> dict[str, Any]:
    """
    Parse {prefix}.mosdepth.summary.txt for global coverage stats.

    The summary file has a header line then rows:
      chrom  length  bases  mean  min  max
    The last row for 'total' gives whole-genome stats.
    """
    global_mean        = 0.0
    pct_above_20x      = 0.0
    pct_above_10x      = 0.0
    total_region_bases = 0

    if not path.exists():
        # Fall back to computing from per_region
        if per_region:
            depths = [r["mean_depth"] for r in per_region]
            global_mean = sum(depths) / len(depths)
        return {
            "mean_coverage":      round(global_mean, 2),
            "pct_above_20x":      0.0,
            "pct_above_10x":      0.0,
            "total_region_bases": 0,
        }

    with open(path) as fh:
        header = None
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if header is None:
                header = parts
                continue
            # chrom  length  bases  mean  min  max
            # "total_region" row has length = total BED target bases
            if parts[0].lower() == "total_region" and len(parts) > 3:
                global_mean        = safe_float(parts[3])
                total_region_bases = safe_int(parts[1])
            elif parts[0].lower() == "total" and global_mean == 0.0 and len(parts) > 3:
                # Fall back to whole-genome total only if no total_region row found yet
                global_mean = safe_float(parts[3])

    # Estimate pct_above thresholds from per_region mean depths.
    # (Exact per-base percentages require per-base parsing — handled in sparklines.)
    # Use a simple heuristic: fraction of regions with mean >= threshold.
    if per_region:
        n = len(per_region)
        n_above_20 = sum(1 for r in per_region if r["mean_depth"] >= min_dp)
        n_above_10 = sum(1 for r in per_region if r["mean_depth"] >= min_dp2)
        pct_above_20x = round(100.0 * n_above_20 / n, 2)
        pct_above_10x = round(100.0 * n_above_10 / n, 2)

        # Backfill pct_20x / pct_10x per region using simple mean-based estimate
        for r in per_region:
            md = r["mean_depth"]
            # Heuristic: Poisson approximation — not exact but good enough without
            # per-base data. We use: P(X >= k) ≈ 1 - Poisson_CDF(k-1, lambda=mean)
            r["pct_20x"] = _estimate_pct_above(md, min_dp)
            r["pct_10x"] = _estimate_pct_above(md, min_dp2)
            r["status"]  = assign_gene_status(md, r["pct_20x"], min_dp)

    return {
        "mean_coverage":      round(global_mean, 2),
        "pct_above_20x":      pct_above_20x,
        "pct_above_10x":      pct_above_10x,
        "total_region_bases": total_region_bases,
    }


def _estimate_pct_above(mean_depth: float, threshold: int) -> float:
    """
    Estimate percent of bases above *threshold* given *mean_depth*.

    Uses a simple sigmoid-like heuristic that gives:
    - 100% when mean >> threshold
    - ~50% when mean == threshold
    - 0%  when mean << threshold
    """
    import math
    if mean_depth <= 0:
        return 0.0
    ratio = mean_depth / threshold
    if ratio >= 3.0:
        return 99.0
    if ratio <= 0.1:
        return 0.0
    # Logistic approximation
    x = (ratio - 1.0) * 3.0
    pct = 100.0 / (1.0 + math.exp(-x))
    return round(min(99.9, max(0.0, pct)), 2)


# ---------------------------------------------------------------------------
# Region aggregation (for BEDs without gene names)
# ---------------------------------------------------------------------------

def _aggregate_regions(
    per_region: list[dict[str, Any]],
    min_dp: int,
    min_dp2: int = 10,
) -> list[dict[str, Any]]:
    """
    Aggregate per-interval records when the BED has no gene names.

    A 3-column BED gives coordinate-style gene labels (e.g. "chr1:65564-65573").
    With 287 K intervals the report becomes 55 MB and unusable.

    Strategy:
    - If ALL gene labels look like "chrom:start-end" (coordinate-style), group
      intervals by chromosome and compute a length-weighted mean depth.
    - Otherwise (real gene names present) return unchanged.

    The result has one record per chromosome (≤ 25 records for a human exome),
    small enough to render instantly.
    """
    if not per_region:
        return per_region

    # Detect whether we have real gene names
    coord_pat_count = sum(
        1 for r in per_region
        if ":" in r.get("gene", "") and "-" in r.get("gene", "")
    )
    has_real_genes = coord_pat_count < len(per_region) * 0.9  # <90% coord-style

    if has_real_genes:
        return per_region  # already has real gene names → return as-is

    # Group by chromosome — length-weighted mean depth
    from collections import defaultdict
    chrom_bases:  dict[str, int]   = defaultdict(int)
    chrom_depth:  dict[str, float] = defaultdict(float)
    chrom_start:  dict[str, int]   = {}
    chrom_end:    dict[str, int]   = {}

    for r in per_region:
        chrom  = r["chrom"]
        length = max(1, r["end"] - r["start"])
        chrom_bases[chrom]  += length
        chrom_depth[chrom]  += r["mean_depth"] * length
        if chrom not in chrom_start:
            chrom_start[chrom] = r["start"]
        chrom_end[chrom] = r["end"]

    aggregated: list[dict[str, Any]] = []
    for chrom in sorted(chrom_bases.keys(), key=_chrom_sort_key):
        total_bases  = chrom_bases[chrom]
        mean_depth   = round(chrom_depth[chrom] / total_bases, 2) if total_bases > 0 else 0.0
        pct_20x      = _estimate_pct_above(mean_depth, min_dp)
        pct_10x      = _estimate_pct_above(mean_depth, min_dp2)
        n_intervals  = len([r for r in per_region if r["chrom"] == chrom])
        start        = chrom_start[chrom]
        end          = chrom_end[chrom]
        # Single-interval chromosomes keep the coordinate label for clarity;
        # multi-interval chromosomes are labelled by chromosome name.
        gene_label   = f"{chrom}:{start}-{end}" if n_intervals == 1 else chrom

        aggregated.append({
            "chrom":       chrom,
            "start":       start,
            "end":         end,
            "gene":        gene_label,
            "mean_depth":  mean_depth,
            "pct_20x":     pct_20x,
            "pct_10x":     pct_10x,
            "status":      assign_gene_status(mean_depth, pct_20x, min_dp),
            "n_intervals": n_intervals,
        })

    return aggregated


def _chrom_sort_key(chrom: str):
    """Natural sort key for chromosome names (chr1 < chr2 … < chrX < chrY < chrM)."""
    name = chrom.replace("chr", "").replace("chrom", "")
    order = {"X": 23, "Y": 24, "M": 25, "MT": 25}
    if name in order:
        return (order[name], 0)
    try:
        return (int(name), 0)
    except ValueError:
        return (99, 0)


# ---------------------------------------------------------------------------
# SIF binary discovery
# ---------------------------------------------------------------------------

def _find_mosdepth_in_sif(singularity_exe: str, sif_path: str) -> str:
    """
    Return the path to the mosdepth binary inside *sif_path*.

    Checks common locations (/opt/mosdepth, /usr/local/bin/mosdepth, mosdepth)
    and returns the first one that exists, falling back to 'mosdepth' (PATH lookup).
    """
    candidates = ["/opt/mosdepth", "/usr/local/bin/mosdepth", "/usr/bin/mosdepth"]
    for candidate in candidates:
        try:
            result = subprocess.run(
                [singularity_exe, "exec", sif_path, "test", "-x", candidate],
                capture_output=True,
            )
            if result.returncode == 0:
                return candidate
        except Exception:
            pass
    # Fall back to PATH lookup inside container
    return "mosdepth"


# ---------------------------------------------------------------------------
# Per-base sparkline index
# ---------------------------------------------------------------------------

def _build_sparkline_index(perbase_file: Path) -> dict[str, list]:
    """
    Read per-base BED and build a look-up dict:
      "{chrom}:{pos}" -> list of (pos, depth) for a ±100bp window

    We store ALL positions and let the JS trim to the 200bp window.
    For very large files we only keep positions that fall inside the
    regions BED (already limited by mosdepth --by).
    """
    index: dict[str, list] = {}  # key: "chrom:start-end" region label

    # We accumulate per-base depths indexed by chrom and position.
    # Then report.py serialises only the 200bp window around each variant.
    # Here we just return the raw stream as a nested dict for slicing later.
    chrom_depths: dict[str, dict[int, int]] = {}

    with gzip.open(perbase_file, "rt") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            chrom = parts[0]
            pos   = safe_int(parts[1])
            depth = safe_int(parts[3])
            if chrom not in chrom_depths:
                chrom_depths[chrom] = {}
            chrom_depths[chrom][pos] = depth

    return chrom_depths  # type: ignore[return-value]


# ---------------------------------------------------------------------------
# Sparkline window extraction
# ---------------------------------------------------------------------------

def get_sparkline_window(
    sparkline_data: dict[str, Any],
    chrom: str,
    pos: int,
    window: int = 100,
) -> list[dict[str, int]]:
    """
    Extract per-base depth for [pos-window, pos+window] from sparkline_data.

    Returns list of {"pos": int, "depth": int} sorted by position.
    """
    chrom_data = sparkline_data.get(chrom, {})
    start = max(0, pos - window)
    end   = pos + window

    result = []
    for p in range(start, end + 1):
        depth = chrom_data.get(p, 0)
        result.append({"pos": p, "depth": depth})
    return result
