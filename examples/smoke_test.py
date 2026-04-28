#!/usr/bin/env python3
"""
smoke_test.py — End-to-end smoke test for vcfdash.

Usage (from repo root):
    python examples/smoke_test.py \
        --bam  /path/to/sample.sorted.bam \
        --vcf  /path/to/sample.vep.vcf.gz \
        --bed  /path/to/panel.bed \
        --sample-id SAMPLE01 \
        --outdir /path/to/output

Requires mosdepth in PATH or --mosdepth-sif <path>.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from vcfdash.coverage import parse_mosdepth_output, run_mosdepth
from vcfdash.variants import parse_vcf
from vcfdash.report import compute_qc, render_html, export_json


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam",         required=True, help="Sorted, indexed BAM")
    parser.add_argument("--vcf",         required=True, help="Annotated VCF (.vcf.gz)")
    parser.add_argument("--bed",         required=True, help="Target BED file")
    parser.add_argument("--sample-id",   required=True, dest="sample_id", help="Sample identifier")
    parser.add_argument("--outdir",      default=".", help="Output directory (default: .)")
    parser.add_argument("--genome",      default="hg38", choices=["hg38", "hg19", "GRCh38", "GRCh37"])
    parser.add_argument("--min-dp",      type=int, default=20, dest="min_dp")
    parser.add_argument("--min-dp2",     type=int, default=10, dest="min_dp2")
    parser.add_argument("--json",        action="store_true", default=False)
    parser.add_argument("--no-sparklines", action="store_true", default=False, dest="no_sparklines")
    parser.add_argument("--mosdepth-sif", default=None, dest="mosdepth_sif",
                        help="Path to mosdepth Singularity SIF image")
    args = parser.parse_args()

    print(f"\nvcfdash smoke test — {args.sample_id}")
    print("=" * 60)

    # ── Step 1: Run mosdepth ─────────────────────────────────────
    print("Running mosdepth …")
    coverage_data = run_mosdepth(
        bam           = args.bam,
        bed           = args.bed,
        outdir        = args.outdir,
        prefix        = args.sample_id,
        min_dp        = args.min_dp,
        min_dp2       = args.min_dp2,
        no_sparklines = args.no_sparklines,
        mosdepth_sif  = args.mosdepth_sif,
    )

    # ── Step 2: Parse VCF ─────────────────────────────────────────
    print("Parsing VCF …")
    variant_data, vcf_warnings = parse_vcf(args.vcf, no_sparklines=args.no_sparklines)

    # ── Step 3: Compute QC ────────────────────────────────────────
    qc = compute_qc(coverage_data, variant_data, args)

    # ── Step 4: Render HTML ───────────────────────────────────────
    html_path = render_html(
        qc_summary    = qc,
        coverage_data = coverage_data,
        variant_data  = variant_data,
        args          = args,
        warnings      = vcf_warnings,
    )

    # ── Step 5: JSON export (optional) ───────────────────────────
    json_path = None
    if args.json:
        json_path = export_json(qc, coverage_data, variant_data, args)

    # ── Summary ───────────────────────────────────────────────────
    per_region = coverage_data.get("per_region", [])
    n_pass = sum(1 for r in per_region if r.get("status") == "PASS")
    n_warn = sum(1 for r in per_region if r.get("status") == "WARN")
    n_fail = sum(1 for r in per_region if r.get("status") == "FAIL")

    import os
    html_mb = os.path.getsize(html_path) / 1_048_576

    print(f"\n{'─'*60}")
    print(f"  vcfdash report complete — {args.sample_id}")
    print(f"{'─'*60}")
    print(f"  Mean coverage (on-BED): {qc.get('mean_coverage', 0):.1f}x  "
          f"({qc.get('bed_target_mb', 0)} Mb BED)")
    print(f"  Bases ≥{args.min_dp}x:          {qc.get('pct_bases_above_20x', 0):.1f}%")
    print(f"  Total variants:         {qc.get('total_variants', 0)}"
          f"  (SNV {qc.get('snv_count',0)}, INDEL {qc.get('indel_count',0)})")
    print(f"  Ti/Tv ratio:            {qc.get('titv_ratio', 0):.2f}")
    print(f"  Regions:                {len(per_region)}"
          f"  (PASS {n_pass}, WARN {n_warn}, FAIL {n_fail})")
    print(f"{'─'*60}")
    print(f"  HTML report:  {html_path}  ({html_mb:.1f} MB)")
    if json_path:
        json_mb = os.path.getsize(json_path) / 1_048_576
        print(f"  JSON summary: {json_path}  ({json_mb:.1f} MB)")
    print(f"{'─'*60}\n")


if __name__ == "__main__":
    main()
