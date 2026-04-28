"""
test_report.py — Tests for vcfdash.report module.
"""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from vcfdash.report import (
    compute_qc,
    export_json,
    embed_assets,
    _build_js_payload,
    _output_path,
)


# ---------------------------------------------------------------------------
# compute_qc
# ---------------------------------------------------------------------------

class TestComputeQc:
    def test_total_variants(self, coverage_data, variant_data, args):
        qc = compute_qc(coverage_data, variant_data, args)
        assert qc["total_variants"] == len(variant_data)

    def test_snv_count(self, coverage_data, variant_data, args):
        qc = compute_qc(coverage_data, variant_data, args)
        expected = sum(1 for v in variant_data if v["variant_type"] == "SNV")
        assert qc["snv_count"] == expected

    def test_indel_count(self, coverage_data, variant_data, args):
        qc = compute_qc(coverage_data, variant_data, args)
        expected = sum(1 for v in variant_data if v["variant_type"] == "INDEL")
        assert qc["indel_count"] == expected

    def test_pass_count(self, coverage_data, variant_data, args):
        qc = compute_qc(coverage_data, variant_data, args)
        expected = sum(1 for v in variant_data if v.get("filter", "PASS") in ("PASS", ".", ""))
        assert qc["pass_variants"] == expected

    def test_flagged_count(self, coverage_data, variant_data, args):
        qc = compute_qc(coverage_data, variant_data, args)
        assert qc["flagged_variants"] == qc["total_variants"] - qc["pass_variants"]

    def test_titv_ratio_positive(self, coverage_data, variant_data, args):
        qc = compute_qc(coverage_data, variant_data, args)
        assert qc["titv_ratio"] >= 0.0

    def test_mean_coverage_from_global_stats(self, coverage_data, args):
        qc = compute_qc(coverage_data, [], args)
        # global stats mean is 59.33 in our synthetic data
        assert qc["mean_coverage"] == pytest.approx(59.33, abs=0.1)

    def test_all_keys_present(self, coverage_data, variant_data, args):
        qc = compute_qc(coverage_data, variant_data, args)
        required = [
            "mean_coverage", "pct_bases_above_20x", "pct_bases_above_10x",
            "total_variants", "snv_count", "indel_count", "mnv_count",
            "titv_ratio", "pass_variants", "flagged_variants",
        ]
        for k in required:
            assert k in qc, f"Missing QC key: {k}"

    def test_empty_variants(self, coverage_data, args):
        qc = compute_qc(coverage_data, [], args)
        assert qc["total_variants"] == 0
        assert qc["snv_count"] == 0
        assert qc["titv_ratio"] == 0.0


# ---------------------------------------------------------------------------
# export_json
# ---------------------------------------------------------------------------

class TestExportJson:
    def test_json_file_created(self, coverage_data, variant_data, args):
        from pathlib import Path
        Path(args.outdir).mkdir(parents=True, exist_ok=True)
        qc = compute_qc(coverage_data, variant_data, args)
        path = export_json(qc, coverage_data, variant_data, args)
        assert Path(path).exists()

    def test_json_schema_keys(self, coverage_data, variant_data, args):
        from pathlib import Path
        Path(args.outdir).mkdir(parents=True, exist_ok=True)
        qc = compute_qc(coverage_data, variant_data, args)
        path = export_json(qc, coverage_data, variant_data, args)
        data = json.loads(Path(path).read_text())
        for key in ("sample_id", "genome", "generated_at", "thresholds",
                    "qc_summary", "gene_coverage", "variants"):
            assert key in data, f"Missing JSON key: {key}"

    def test_json_sample_id(self, coverage_data, variant_data, args):
        from pathlib import Path
        Path(args.outdir).mkdir(parents=True, exist_ok=True)
        qc = compute_qc(coverage_data, variant_data, args)
        path = export_json(qc, coverage_data, variant_data, args)
        data = json.loads(Path(path).read_text())
        assert data["sample_id"] == "TEST_SAMPLE"

    def test_json_variants_count(self, coverage_data, variant_data, args):
        from pathlib import Path
        Path(args.outdir).mkdir(parents=True, exist_ok=True)
        qc = compute_qc(coverage_data, variant_data, args)
        path = export_json(qc, coverage_data, variant_data, args)
        data = json.loads(Path(path).read_text())
        assert len(data["variants"]) == len(variant_data)

    def test_json_gene_coverage_count(self, coverage_data, variant_data, args):
        from pathlib import Path
        Path(args.outdir).mkdir(parents=True, exist_ok=True)
        qc = compute_qc(coverage_data, variant_data, args)
        path = export_json(qc, coverage_data, variant_data, args)
        data = json.loads(Path(path).read_text())
        expected = len(coverage_data.get("per_region", []))
        assert len(data["gene_coverage"]) == expected

    def test_json_gene_coverage_status_field(self, coverage_data, variant_data, args):
        from pathlib import Path
        Path(args.outdir).mkdir(parents=True, exist_ok=True)
        qc = compute_qc(coverage_data, variant_data, args)
        path = export_json(qc, coverage_data, variant_data, args)
        data = json.loads(Path(path).read_text())
        for gc in data["gene_coverage"]:
            assert gc["status"] in ("PASS", "WARN", "FAIL")

    def test_json_variant_fields(self, coverage_data, variant_data, args):
        from pathlib import Path
        Path(args.outdir).mkdir(parents=True, exist_ok=True)
        qc = compute_qc(coverage_data, variant_data, args)
        path = export_json(qc, coverage_data, variant_data, args)
        data = json.loads(Path(path).read_text())
        v = data["variants"][0]
        for field in ("chrom", "pos", "ref", "alt", "gene", "consequence",
                      "hgvs_c", "hgvs_p", "gnomad_af", "clinvar",
                      "cadd_phred", "gt", "dp", "gq", "vaf", "filter"):
            assert field in v, f"Missing variant field in JSON: {field}"

    def test_json_thresholds(self, coverage_data, variant_data, args):
        from pathlib import Path
        Path(args.outdir).mkdir(parents=True, exist_ok=True)
        qc = compute_qc(coverage_data, variant_data, args)
        path = export_json(qc, coverage_data, variant_data, args)
        data = json.loads(Path(path).read_text())
        assert data["thresholds"]["min_dp"] == 20
        assert data["thresholds"]["min_dp2"] == 10

    def test_json_is_valid(self, coverage_data, variant_data, args):
        """The output must be valid JSON."""
        from pathlib import Path
        Path(args.outdir).mkdir(parents=True, exist_ok=True)
        qc = compute_qc(coverage_data, variant_data, args)
        path = export_json(qc, coverage_data, variant_data, args)
        raw = Path(path).read_text()
        parsed = json.loads(raw)  # must not raise
        assert isinstance(parsed, dict)


# ---------------------------------------------------------------------------
# embed_assets
# ---------------------------------------------------------------------------

class TestEmbedAssets:
    def test_returns_nonempty_strings(self):
        from vcfdash.report import _TEMPLATE_DIR
        css, viz_js = embed_assets(_TEMPLATE_DIR)
        assert len(css) > 100
        assert len(viz_js) > 100

    def test_css_contains_root_variables(self):
        from vcfdash.report import _TEMPLATE_DIR
        css, _ = embed_assets(_TEMPLATE_DIR)
        assert ":root" in css

    def test_vizjs_contains_vcfdash_key(self):
        from vcfdash.report import _TEMPLATE_DIR
        _, viz_js = embed_assets(_TEMPLATE_DIR)
        assert "VCFDASH" in viz_js


# ---------------------------------------------------------------------------
# _build_js_payload
# ---------------------------------------------------------------------------

class TestBuildJsPayload:
    def test_returns_valid_json(self, coverage_data, variant_data, args):
        payload_str = _build_js_payload({}, coverage_data, variant_data, args)
        parsed = json.loads(payload_str)
        assert isinstance(parsed, dict)

    def test_contains_required_keys(self, coverage_data, variant_data, args):
        payload_str = _build_js_payload({}, coverage_data, variant_data, args)
        parsed = json.loads(payload_str)
        # New sparkline format: parallel arrays, not a full genome dict
        for k in ("sample_id", "thresholds", "qc", "coverage", "variants",
                  "spark_windows", "spark_start", "has_sparklines"):
            assert k in parsed, f"Missing payload key: {k}"

    def test_variants_present(self, coverage_data, variant_data, args):
        payload_str = _build_js_payload({}, coverage_data, variant_data, args)
        parsed = json.loads(payload_str)
        assert len(parsed["variants"]) == len(variant_data)

    def test_sparklines_empty_when_no_sparklines(self, coverage_data, variant_data, args):
        args.no_sparklines = True
        payload_str = _build_js_payload({}, coverage_data, variant_data, args)
        parsed = json.loads(payload_str)
        # With no_sparklines=True, windows arrays should be empty and flag False
        assert parsed["spark_windows"] == []
        assert parsed["spark_start"]   == []
        assert parsed["has_sparklines"] is False

    def test_ad_serialised_as_list(self, coverage_data, variant_data, args):
        payload_str = _build_js_payload({}, coverage_data, variant_data, args)
        parsed = json.loads(payload_str)
        for v in parsed["variants"]:
            if v.get("ad") is not None:
                assert isinstance(v["ad"], list)


# ---------------------------------------------------------------------------
# render_html (integration — no D3 CDN required)
# ---------------------------------------------------------------------------

class TestRenderHtml:
    def test_html_file_created(self, coverage_data, variant_data, args, monkeypatch):
        """render_html should write a .html file without needing mosdepth or real VCF."""
        from pathlib import Path
        from vcfdash.report import render_html, _D3_FALLBACK
        import vcfdash.report as rmod

        # Stub out D3 fetch so test works offline
        monkeypatch.setattr(rmod, "_D3_CACHE", _D3_FALLBACK)

        Path(args.outdir).mkdir(parents=True, exist_ok=True)
        qc = compute_qc(coverage_data, variant_data, args)
        path = render_html(qc, coverage_data, variant_data, args)
        assert Path(path).exists()
        assert path.endswith(".html")

    def test_html_contains_sample_id(self, coverage_data, variant_data, args, monkeypatch):
        from pathlib import Path
        from vcfdash.report import render_html, _D3_FALLBACK
        import vcfdash.report as rmod

        monkeypatch.setattr(rmod, "_D3_CACHE", _D3_FALLBACK)
        Path(args.outdir).mkdir(parents=True, exist_ok=True)
        qc = compute_qc(coverage_data, variant_data, args)
        path = render_html(qc, coverage_data, variant_data, args)
        html = Path(path).read_text()
        assert "TEST_SAMPLE" in html

    def test_html_is_self_contained(self, coverage_data, variant_data, args, monkeypatch):
        """Output must not contain http:// CDN references (fully self-contained)."""
        from pathlib import Path
        from vcfdash.report import render_html, _D3_FALLBACK
        import vcfdash.report as rmod

        monkeypatch.setattr(rmod, "_D3_CACHE", _D3_FALLBACK)
        Path(args.outdir).mkdir(parents=True, exist_ok=True)
        qc = compute_qc(coverage_data, variant_data, args)
        path = render_html(qc, coverage_data, variant_data, args)
        html = Path(path).read_text()
        # No external links to CDN resources
        assert 'src="http' not in html
        assert "href=\"http" not in html

    def test_html_contains_vcfdash_json_blob(self, coverage_data, variant_data, args, monkeypatch):
        from pathlib import Path
        from vcfdash.report import render_html, _D3_FALLBACK
        import vcfdash.report as rmod

        monkeypatch.setattr(rmod, "_D3_CACHE", _D3_FALLBACK)
        Path(args.outdir).mkdir(parents=True, exist_ok=True)
        qc = compute_qc(coverage_data, variant_data, args)
        path = render_html(qc, coverage_data, variant_data, args)
        html = Path(path).read_text()
        assert "window.VCFDASH" in html

    def test_html_contains_three_panels(self, coverage_data, variant_data, args, monkeypatch):
        from pathlib import Path
        from vcfdash.report import render_html, _D3_FALLBACK
        import vcfdash.report as rmod

        monkeypatch.setattr(rmod, "_D3_CACHE", _D3_FALLBACK)
        Path(args.outdir).mkdir(parents=True, exist_ok=True)
        qc = compute_qc(coverage_data, variant_data, args)
        path = render_html(qc, coverage_data, variant_data, args)
        html = Path(path).read_text()
        # Check for panel markers
        assert "QC Summary" in html
        assert "Coverage" in html
        assert "Variants" in html


# ---------------------------------------------------------------------------
# _output_path
# ---------------------------------------------------------------------------

class TestOutputPath:
    def test_creates_directory(self, args):
        from pathlib import Path
        args.outdir = str(Path(args.outdir) / "subdir")
        p = _output_path(args, "test.html")
        assert p.parent.exists()

    def test_filename_correct(self, args):
        from pathlib import Path
        Path(args.outdir).mkdir(parents=True, exist_ok=True)
        p = _output_path(args, "sample_report.html")
        assert p.name == "sample_report.html"
