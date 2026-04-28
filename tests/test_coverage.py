"""
test_coverage.py — Tests for vcfdash.coverage module.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from vcfdash.coverage import (
    parse_mosdepth_output,
    get_sparkline_window,
    _estimate_pct_above,
)
from vcfdash.utils import assign_gene_status


# ---------------------------------------------------------------------------
# parse_mosdepth_output
# ---------------------------------------------------------------------------

class TestParseMosdepthOutput:
    def test_returns_three_regions(self, mosdepth_prefix):
        data = parse_mosdepth_output(mosdepth_prefix, min_dp=20, min_dp2=10, no_sparklines=True)
        assert len(data["per_region"]) == 3

    def test_region_fields_present(self, mosdepth_prefix):
        data = parse_mosdepth_output(mosdepth_prefix, min_dp=20, min_dp2=10, no_sparklines=True)
        r = data["per_region"][0]
        for key in ("chrom", "start", "end", "gene", "mean_depth", "pct_20x", "pct_10x", "status"):
            assert key in r, f"Missing key: {key}"

    def test_gene_names_extracted(self, mosdepth_prefix):
        data = parse_mosdepth_output(mosdepth_prefix, min_dp=20, min_dp2=10, no_sparklines=True)
        genes = {r["gene"] for r in data["per_region"]}
        assert "GENE_A" in genes
        assert "GENE_B" in genes
        assert "GENE_C" in genes

    def test_mean_depth_values(self, mosdepth_prefix):
        data = parse_mosdepth_output(mosdepth_prefix, min_dp=20, min_dp2=10, no_sparklines=True)
        depths = {r["gene"]: r["mean_depth"] for r in data["per_region"]}
        assert depths["GENE_A"] == pytest.approx(85.3, abs=0.01)
        assert depths["GENE_B"] == pytest.approx(12.1, abs=0.01)
        assert depths["GENE_C"] == pytest.approx(7.4, abs=0.01)

    def test_global_stats_present(self, mosdepth_prefix):
        data = parse_mosdepth_output(mosdepth_prefix, min_dp=20, min_dp2=10, no_sparklines=True)
        gs = data["global_stats"]
        assert "mean_coverage" in gs
        assert "pct_above_20x" in gs
        assert "pct_above_10x" in gs

    def test_global_mean_from_summary(self, mosdepth_prefix):
        data = parse_mosdepth_output(mosdepth_prefix, min_dp=20, min_dp2=10, no_sparklines=True)
        # synthetic summary has total mean of 59.33
        assert data["global_stats"]["mean_coverage"] == pytest.approx(59.33, abs=0.01)

    def test_status_assigned(self, mosdepth_prefix):
        data = parse_mosdepth_output(mosdepth_prefix, min_dp=20, min_dp2=10, no_sparklines=True)
        statuses = {r["gene"]: r["status"] for r in data["per_region"]}
        assert statuses["GENE_A"] == "PASS"     # depth 85.3 >> 20
        assert statuses["GENE_B"] in ("WARN", "FAIL")  # depth 12.1 < 20
        assert statuses["GENE_C"] == "FAIL"     # depth 7.4 << 20

    def test_sparkline_data_empty_when_no_sparklines(self, mosdepth_prefix):
        data = parse_mosdepth_output(mosdepth_prefix, min_dp=20, min_dp2=10, no_sparklines=True)
        assert data["sparkline_data"] == {}

    def test_missing_regions_file_raises(self, tmp_path):
        bad_prefix = str(tmp_path / "nonexistent")
        with pytest.raises(SystemExit):
            parse_mosdepth_output(bad_prefix, min_dp=20, min_dp2=10, no_sparklines=True)

    def test_4col_bed_no_gene_name(self, tmp_path):
        """Regions BED without gene name in col 4 should use chrom:start-end."""
        prefix = str(tmp_path / "sample")
        # Only 4 columns — no gene name
        content = "chr1\t100\t200\t55.0\n"
        with gzip.open(f"{prefix}.regions.bed.gz", "wt") as fh:
            fh.write(content)
        Path(f"{prefix}.mosdepth.summary.txt").write_text(
            "chrom\tlength\tbases\tmean\tmin\tmax\ntotal\t100\t5500\t55.0\t0\t100\n"
        )
        data = parse_mosdepth_output(prefix, min_dp=20, min_dp2=10, no_sparklines=True)
        assert data["per_region"][0]["gene"] == "chr1:100-200"


# ---------------------------------------------------------------------------
# assign_gene_status
# ---------------------------------------------------------------------------

class TestAssignGeneStatus:
    def test_pass_high_depth_high_pct(self):
        assert assign_gene_status(80.0, 97.0, 20) == "PASS"

    def test_warn_medium_depth(self):
        assert assign_gene_status(12.0, 81.0, 20) == "WARN"

    def test_fail_low_depth(self):
        assert assign_gene_status(3.0, 20.0, 20) == "FAIL"

    def test_fail_when_depth_ok_but_pct_low(self):
        # depth >= min_dp but pct_20x < 95 → not PASS
        # depth >= min_dp * 0.5 but pct_20x < 80 → not WARN → FAIL
        assert assign_gene_status(25.0, 50.0, 20) == "FAIL"

    def test_warn_borderline(self):
        # depth = 10 (0.5 * 20), pct = 82 (>= 80) → WARN
        assert assign_gene_status(10.0, 82.0, 20) == "WARN"

    def test_custom_min_dp(self):
        # With min_dp=10, depth=12 >= 10 is PASS threshold
        assert assign_gene_status(12.0, 96.0, 10) == "PASS"


# ---------------------------------------------------------------------------
# _estimate_pct_above
# ---------------------------------------------------------------------------

class TestEstimatePctAbove:
    def test_zero_depth_returns_zero(self):
        assert _estimate_pct_above(0.0, 20) == 0.0

    def test_high_ratio_returns_near_100(self):
        pct = _estimate_pct_above(100.0, 20)  # ratio = 5
        assert pct >= 98.0

    def test_ratio_one_returns_near_50(self):
        pct = _estimate_pct_above(20.0, 20)   # ratio = 1
        assert 40.0 <= pct <= 60.0

    def test_low_ratio_returns_near_zero(self):
        pct = _estimate_pct_above(1.0, 20)    # ratio = 0.05
        assert pct == 0.0


# ---------------------------------------------------------------------------
# get_sparkline_window
# ---------------------------------------------------------------------------

class TestGetSparklineWindow:
    def test_returns_correct_range(self):
        data = {"chr1": {100: 30, 150: 45, 200: 20}}
        window = get_sparkline_window(data, "chr1", 150, window=10)
        positions = [pt["pos"] for pt in window]
        assert min(positions) == 140
        assert max(positions) == 160

    def test_missing_positions_return_zero(self):
        data = {"chr1": {150: 30}}
        window = get_sparkline_window(data, "chr1", 150, window=5)
        for pt in window:
            if pt["pos"] != 150:
                assert pt["depth"] == 0

    def test_known_position_depth(self):
        data = {"chr1": {150: 42}}
        window = get_sparkline_window(data, "chr1", 150, window=5)
        hit = next(pt for pt in window if pt["pos"] == 150)
        assert hit["depth"] == 42

    def test_missing_chrom_returns_zeros(self):
        data = {}
        window = get_sparkline_window(data, "chrX", 1000, window=5)
        assert all(pt["depth"] == 0 for pt in window)

    def test_window_size_respected(self):
        data = {"chr1": {p: 10 for p in range(1, 300)}}
        window = get_sparkline_window(data, "chr1", 150, window=50)
        assert len(window) == 101  # 50 + 1 + 50
