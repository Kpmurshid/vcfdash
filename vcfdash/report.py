"""
report.py — Jinja2 HTML rendering and JSON export.

All CSS, JS, and data are inlined into the output HTML file so it is
fully self-contained and opens without any internet connection or server.
"""

from __future__ import annotations

import json
import sys
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


# ---------------------------------------------------------------------------
# JSON encoder that handles numpy scalar types (int32, float32, etc.)
# ---------------------------------------------------------------------------

class _NumpySafeEncoder(json.JSONEncoder):
    """JSON encoder that converts numpy scalars to native Python types."""
    def default(self, obj):
        # Handle numpy integer types
        try:
            import numpy as np
            if isinstance(obj, np.integer):
                return int(obj)
            if isinstance(obj, np.floating):
                return float(obj)
            if isinstance(obj, np.ndarray):
                return obj.tolist()
            if isinstance(obj, np.bool_):
                return bool(obj)
        except ImportError:
            pass
        return super().default(obj)

from . import __version__
from .coverage import get_sparkline_window
from .variants import compute_titv

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

_TEMPLATE_DIR = Path(__file__).parent / "templates"
_ASSETS_DIR   = _TEMPLATE_DIR / "assets"

# D3.js v7 CDN URL — fetched once at build time and embedded
_D3_CDN_URL   = "https://cdn.jsdelivr.net/npm/d3@7/dist/d3.min.js"
_D3_FALLBACK  = "// D3.js not available — sparklines disabled"

# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def render_html(
    qc_summary: dict,
    coverage_data: dict,
    variant_data: list[dict],
    args,
    warnings: list[str] | None = None,
) -> str:
    """
    Render the HTML report and write it to {outdir}/{sample_id}_report.html.

    Returns the output file path as a string.
    """
    from jinja2 import Environment, FileSystemLoader, select_autoescape

    css, viz_js = embed_assets(_TEMPLATE_DIR)
    d3_js       = _get_d3_js()

    vcfdash_payload = _build_js_payload(
        qc_summary, coverage_data, variant_data, args
    )

    env = Environment(
        loader=FileSystemLoader(str(_TEMPLATE_DIR)),
        autoescape=select_autoescape(["html"]),
        # Disable autoescaping for our known-safe JSON/CSS/JS blobs
    )
    # We control the JSON serialisation ourselves, so disable HTML escaping
    # for the script-tag content by using |safe in the template (which we do).
    env.policies["json.dumps_kwargs"] = {"ensure_ascii": False}

    template = env.get_template("report.html.j2")

    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

    html = template.render(
        sample_id    = args.sample_id,
        genome       = getattr(args, "genome", "hg38"),
        generated_at = now,
        version      = __version__,
        warnings     = warnings or [],
        qc           = _DotDict(qc_summary),
        thresholds   = _DotDict({
            "min_dp":  args.min_dp,
            "min_dp2": args.min_dp2,
        }),
        css          = css,
        d3_js        = d3_js,
        viz_js       = viz_js,
        vcfdash_json = vcfdash_payload,
    )

    out_path = _output_path(args, f"{args.sample_id}_report.html")
    out_path.write_text(html, encoding="utf-8")
    print(f"[vcfdash] HTML report written → {out_path}", file=sys.stderr)
    return str(out_path)


def export_json(
    qc_summary: dict,
    coverage_data: dict,
    variant_data: list[dict],
    args,
) -> str:
    """
    Write {sample_id}_report.json following the specified schema.
    Only called when --json flag is passed.

    Returns the output file path as a string.
    """
    now = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    gene_coverage = [
        {
            "gene":         r.get("gene", ""),
            "mean_depth":   round(r.get("mean_depth", 0.0), 2),
            "pct_above_20x": round(r.get("pct_20x") or 0.0, 2),
            "pct_above_10x": round(r.get("pct_10x") or 0.0, 2),
            "status":        r.get("status", "FAIL"),
        }
        for r in coverage_data.get("per_region", [])
    ]

    variants_export = [
        {
            "chrom":      v.get("chrom", ""),
            "pos":        v.get("pos", 0),
            "ref":        v.get("ref", ""),
            "alt":        v.get("alt", ""),
            "gene":       v.get("gene", ""),
            "consequence": v.get("consequence", ""),
            "hgvs_c":     v.get("hgvs_c", ""),
            "hgvs_p":     v.get("hgvs_p", ""),
            "gnomad_af":  v.get("gnomad_af") or 0.0,
            "clinvar":    v.get("clinvar", ""),
            "cadd_phred": v.get("cadd_phred") or 0.0,
            "gt":         v.get("gt", "."),
            "dp":         v.get("dp") or 0,
            "gq":         v.get("gq") or 0,
            "vaf":        v.get("vaf"),
            "filter":     v.get("filter", "PASS"),
        }
        for v in variant_data
    ]

    payload = {
        "sample_id":    args.sample_id,
        "genome":       getattr(args, "genome", "hg38"),
        "generated_at": now,
        "thresholds": {
            "min_dp":  args.min_dp,
            "min_dp2": args.min_dp2,
        },
        "qc_summary":    qc_summary,
        "gene_coverage": gene_coverage,
        "variants":      variants_export,
    }

    out_path = _output_path(args, f"{args.sample_id}_report.json")
    out_path.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False, cls=_NumpySafeEncoder),
        encoding="utf-8",
    )
    print(f"[vcfdash] JSON summary written  → {out_path}", file=sys.stderr)
    return str(out_path)


def embed_assets(template_dir: Path) -> tuple[str, str]:
    """
    Read CSS and JS assets from disk and return them as strings for inline embedding.

    Returns (css_string, viz_js_string).
    """
    css_path    = template_dir / "assets" / "style.css"
    viz_js_path = template_dir / "assets" / "viz.js"

    css    = css_path.read_text(encoding="utf-8") if css_path.exists() else ""
    viz_js = viz_js_path.read_text(encoding="utf-8") if viz_js_path.exists() else ""

    return css, viz_js


# ---------------------------------------------------------------------------
# JS payload builder
# ---------------------------------------------------------------------------

def _build_js_payload(
    qc_summary: dict,
    coverage_data: dict,
    variant_data: list[dict],
    args,
) -> str:
    """
    Serialise all report data to a JSON string for inline embedding.

    Sparkline data is included only when available and --no-sparklines not set.
    """
    no_sparklines  = getattr(args, "no_sparklines", False)
    sparkline_data = coverage_data.get("sparkline_data", {})

    # Build per-variant sparkline windows
    # We serialize chrom -> {pos: depth} for the whole covered region and let
    # viz.js slice the 200bp window on the client side.
    # For large datasets this may be megabytes; use no_sparklines to disable.
    spark_payload: dict[str, Any] = {}
    if not no_sparklines and sparkline_data:
        # Convert position keys to strings for JSON compatibility
        for chrom, positions in sparkline_data.items():
            spark_payload[chrom] = {str(p): d for p, d in positions.items()}

    # Sanitize variants for JSON (remove non-serializable ad arrays → lists)
    variants_clean = []
    for v in variant_data:
        vc = dict(v)
        if isinstance(vc.get("ad"), (list, tuple)):
            vc["ad"] = list(vc["ad"])
        variants_clean.append(vc)

    payload = {
        "sample_id":      args.sample_id,
        "genome":         getattr(args, "genome", "hg38"),
        "has_sparklines": not no_sparklines and bool(sparkline_data),
        "thresholds": {
            "min_dp":  args.min_dp,
            "min_dp2": args.min_dp2,
        },
        "qc":        qc_summary,
        "coverage":  coverage_data.get("per_region", []),
        "variants":  variants_clean,
        "sparklines": spark_payload,
    }

    return json.dumps(payload, ensure_ascii=False, allow_nan=False, cls=_NumpySafeEncoder)


# ---------------------------------------------------------------------------
# D3.js fetching / caching
# ---------------------------------------------------------------------------

_D3_CACHE: str | None = None


def _get_d3_js() -> str:
    """
    Return the D3.js v7 source code string.

    Tries (in order):
    1. Local cache (same Python session)
    2. Bundled fallback in assets/d3.min.js  (if pre-downloaded)
    3. Download from CDN
    4. Return a stub that disables sparklines gracefully
    """
    global _D3_CACHE
    if _D3_CACHE is not None:
        return _D3_CACHE

    # Check for pre-bundled copy
    bundled = _ASSETS_DIR / "d3.min.js"
    if bundled.exists():
        _D3_CACHE = bundled.read_text(encoding="utf-8")
        return _D3_CACHE

    # Try to download
    print("[vcfdash] Fetching D3.js v7 from CDN …", file=sys.stderr)
    try:
        with urllib.request.urlopen(_D3_CDN_URL, timeout=15) as resp:
            _D3_CACHE = resp.read().decode("utf-8")
            # Save for next time
            try:
                bundled.write_text(_D3_CACHE, encoding="utf-8")
            except OSError:
                pass
            return _D3_CACHE
    except Exception as exc:
        print(
            f"[vcfdash] WARNING: Could not fetch D3.js ({exc}). "
            "Sparklines will be disabled in output HTML.",
            file=sys.stderr,
        )
        _D3_CACHE = _D3_FALLBACK
        return _D3_CACHE


# ---------------------------------------------------------------------------
# QC summary builder
# ---------------------------------------------------------------------------

def compute_qc(
    coverage_data: dict,
    variant_data: list[dict],
    args,
) -> dict[str, Any]:
    """
    Compute the QC summary dict from coverage + variant data.
    """
    global_stats = coverage_data.get("global_stats", {})

    total     = len(variant_data)
    snv_count = sum(1 for v in variant_data if v.get("variant_type") == "SNV")
    ind_count = sum(1 for v in variant_data if v.get("variant_type") == "INDEL")
    mnv_count = sum(1 for v in variant_data if v.get("variant_type") == "MNV")
    pass_count = sum(
        1 for v in variant_data
        if (v.get("filter") or "PASS") in ("PASS", ".", "")
    )
    flagged   = total - pass_count
    titv      = compute_titv(variant_data)

    # BED target size — read from the mosdepth summary's total_region row
    # (the "length" column = sum of all interval lengths in the BED file).
    # This is the authoritative count of targeted bases and matches what
    # Picard HsMetrics calls "TARGET_TERRITORY".
    bed_bases     = global_stats.get("total_region_bases", 0)
    bed_target_mb = round(bed_bases / 1_000_000, 1) if bed_bases else 0.0

    return {
        "mean_coverage":       global_stats.get("mean_coverage", 0.0),
        "pct_bases_above_20x": global_stats.get("pct_above_20x", 0.0),
        "pct_bases_above_10x": global_stats.get("pct_above_10x", 0.0),
        "total_variants":      total,
        "snv_count":           snv_count,
        "indel_count":         ind_count,
        "mnv_count":           mnv_count,
        "titv_ratio":          titv,
        "pass_variants":       pass_count,
        "flagged_variants":    flagged,
        "bed_target_mb":       bed_target_mb,
    }


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _output_path(args, filename: str) -> Path:
    outdir = Path(getattr(args, "outdir", ".") or ".")
    outdir.mkdir(parents=True, exist_ok=True)
    return outdir / filename


class _DotDict(dict):
    """Dict subclass that allows attribute-style access in Jinja2 templates."""
    def __getattr__(self, name):
        try:
            return self[name]
        except KeyError:
            raise AttributeError(name)
