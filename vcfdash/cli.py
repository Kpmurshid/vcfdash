"""
cli.py — Command-line entry point for vcfdash.

Orchestrates the pipeline:
  1. Validate inputs
  2. Run mosdepth → parse coverage
  3. Parse VCF → extract variants
  4. Compute QC summary
  5. Render HTML report
  6. Optionally write JSON export
  7. Print summary to stdout
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import __version__


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="vcfdash",
        description=(
            "vcfdash — Lightweight, pipeline-native clinical variant and coverage dashboard.\n"
            "Produces a self-contained HTML report from BAM + VCF + BED inputs."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  vcfdash --bam sample.bam --vcf sample.vep.vcf.gz --bed panel.bed --sample-id SAMPLE01
  vcfdash --bam sample.bam --vcf sample.vcf --bed panel.bed --sample-id S1 --json --outdir reports/
        """,
    )

    parser.add_argument("--version", action="version", version=f"vcfdash {__version__}")

    # Required
    req = parser.add_argument_group("required arguments")
    req.add_argument(
        "--bam", required=True, metavar="PATH",
        help="Sorted, indexed BAM file (.bam + .bai must exist)",
    )
    req.add_argument(
        "--vcf", required=True, metavar="PATH",
        help="Annotated VCF file (post-VEP or post-Annovar preferred)",
    )
    req.add_argument(
        "--bed", required=True, metavar="PATH",
        help="Target regions BED file (gene panel or exome capture)",
    )
    req.add_argument(
        "--sample-id", required=True, metavar="STRING", dest="sample_id",
        help="Sample identifier string (used in report title and filenames)",
    )

    # Optional
    opt = parser.add_argument_group("optional arguments")
    opt.add_argument(
        "--outdir", metavar="PATH", default=".",
        help="Output directory (default: current directory)",
    )
    opt.add_argument(
        "--genome", metavar="STRING", default="hg38",
        choices=["hg38", "hg19", "GRCh38", "GRCh37"],
        help="Reference genome build: hg38 (default) or hg19",
    )
    opt.add_argument(
        "--min-dp", metavar="INT", type=int, default=20, dest="min_dp",
        help="Minimum depth threshold for coverage flagging (default: 20)",
    )
    opt.add_argument(
        "--min-dp2", metavar="INT", type=int, default=10, dest="min_dp2",
        help="Secondary depth threshold (default: 10)",
    )
    opt.add_argument(
        "--config", metavar="PATH", default=None,
        help="JSON config file for custom thresholds and panel metadata",
    )
    opt.add_argument(
        "--json", action="store_true", default=False,
        help="Export machine-readable JSON summary alongside HTML report",
    )
    opt.add_argument(
        "--no-sparklines", action="store_true", default=False, dest="no_sparklines",
        help="Disable per-variant coverage sparklines (faster, smaller output)",
    )
    opt.add_argument(
        "--mosdepth-sif", metavar="PATH", default=None, dest="mosdepth_sif",
        help=(
            "Path to a mosdepth Singularity/Apptainer SIF image. "
            "Use when mosdepth is not installed in PATH. "
            "Auto-detected from well-known locations if omitted."
        ),
    )

    return parser.parse_args(argv)


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)

    # ── Step 1: Validate inputs ────────────────────────────────────────────
    from .utils import (
        check_file_readable,
        check_bam_index,
        check_mosdepth,
        check_output_dir,
        load_config,
    )

    print(f"[vcfdash] vcfdash {__version__} starting …", file=sys.stderr)
    print(f"[vcfdash] Sample: {args.sample_id}", file=sys.stderr)

    check_file_readable(args.bam, "BAM file")
    check_bam_index(args.bam)
    check_file_readable(args.vcf, "VCF file")
    check_file_readable(args.bed, "BED file")
    check_mosdepth(getattr(args, "mosdepth_sif", None))  # validates SIF/PATH early
    check_output_dir(args.outdir)

    # Load config (merges file + CLI args)
    cfg = load_config(args.config, args)

    # ── Step 2: Run mosdepth → parse coverage ────────────────────────────
    from .coverage import run_mosdepth

    coverage_data = run_mosdepth(
        bam           = args.bam,
        bed           = args.bed,
        outdir        = args.outdir,
        prefix        = args.sample_id,
        min_dp        = args.min_dp,
        min_dp2       = args.min_dp2,
        no_sparklines = args.no_sparklines,
        mosdepth_sif  = getattr(args, "mosdepth_sif", None),
    )

    n_regions = len(coverage_data.get("per_region", []))
    print(f"[vcfdash] Coverage parsed: {n_regions} regions", file=sys.stderr)

    # ── Step 3: Parse VCF ────────────────────────────────────────────────
    from .variants import parse_vcf

    variant_data, vcf_warnings = parse_vcf(args.vcf, no_sparklines=args.no_sparklines)
    print(f"[vcfdash] Variants parsed:  {len(variant_data)} records", file=sys.stderr)

    # ── Step 4: Compute QC summary ────────────────────────────────────────
    from .report import compute_qc

    qc_summary = compute_qc(coverage_data, variant_data, args)

    # ── Step 5: Render HTML report ────────────────────────────────────────
    from .report import render_html

    html_path = render_html(
        qc_summary    = qc_summary,
        coverage_data = coverage_data,
        variant_data  = variant_data,
        args          = args,
        warnings      = vcf_warnings,
    )

    # ── Step 6: JSON export (optional) ───────────────────────────────────
    json_path: str | None = None
    if args.json:
        from .report import export_json
        json_path = export_json(qc_summary, coverage_data, variant_data, args)

    # ── Step 7: Print summary to stdout ──────────────────────────────────
    _print_summary(qc_summary, coverage_data, variant_data, html_path, json_path, args)


# ---------------------------------------------------------------------------
# Summary printer
# ---------------------------------------------------------------------------

def _print_summary(
    qc: dict,
    coverage_data: dict,
    variant_data: list,
    html_path: str,
    json_path: str | None,
    args,
) -> None:
    per_region = coverage_data.get("per_region", [])
    n_pass  = sum(1 for r in per_region if r.get("status") == "PASS")
    n_warn  = sum(1 for r in per_region if r.get("status") == "WARN")
    n_fail  = sum(1 for r in per_region if r.get("status") == "FAIL")

    sep = "─" * 60
    print(f"\n{sep}")
    print(f"  vcfdash report complete — {args.sample_id}")
    print(sep)
    print(f"  Mean coverage:     {qc.get('mean_coverage', 0):.1f}x")
    print(f"  Bases ≥{args.min_dp}x:      {qc.get('pct_bases_above_20x', 0):.1f}%")
    print(f"  Bases ≥{args.min_dp2}x:      {qc.get('pct_bases_above_10x', 0):.1f}%")
    print(f"  Total variants:    {qc.get('total_variants', 0)}"
          f"  (SNV {qc.get('snv_count',0)}, INDEL {qc.get('indel_count',0)}, MNV {qc.get('mnv_count',0)})")
    print(f"  Ti/Tv ratio:       {qc.get('titv_ratio', 0):.2f}")
    print(f"  Gene regions:      {len(per_region)}"
          f"  (PASS {n_pass}, WARN {n_warn}, FAIL {n_fail})")
    print(sep)
    print(f"  HTML report:  {html_path}")
    if json_path:
        print(f"  JSON summary: {json_path}")
    print(sep)
    print()


if __name__ == "__main__":
    main()
