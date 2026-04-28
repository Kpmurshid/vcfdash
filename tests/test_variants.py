"""
test_variants.py — Tests for vcfdash.variants module.
"""

from __future__ import annotations

import textwrap
from pathlib import Path

import pytest

from vcfdash.variants import (
    classify_variant_type,
    compute_titv,
    extract_vep_csq,
    _is_more_severe,
    _empty_annotation,
    _parse_vep_format_string,
)


# ---------------------------------------------------------------------------
# classify_variant_type
# ---------------------------------------------------------------------------

class TestClassifyVariantType:
    def test_snv(self):
        assert classify_variant_type("A", "G") == "SNV"
        assert classify_variant_type("C", "T") == "SNV"

    def test_indel_deletion(self):
        assert classify_variant_type("ATG", "A") == "INDEL"

    def test_indel_insertion(self):
        assert classify_variant_type("A", "ATGC") == "INDEL"

    def test_mnv(self):
        assert classify_variant_type("AT", "GC") == "MNV"
        assert classify_variant_type("AAA", "TTT") == "MNV"

    def test_ref_dot(self):
        # spanning deletion
        assert classify_variant_type("A", "*") == "SNV"

    def test_single_base_is_snv(self):
        assert classify_variant_type("G", "A") == "SNV"


# ---------------------------------------------------------------------------
# compute_titv
# ---------------------------------------------------------------------------

class TestComputeTiTv:
    def _make_variants(self, snv_pairs: list[tuple[str, str]]) -> list[dict]:
        return [
            {"ref": r, "alt": a, "variant_type": "SNV"}
            for r, a in snv_pairs
        ]

    def test_all_transitions(self):
        variants = self._make_variants([("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")])
        # 4 Ti, 0 Tv → sentinel 999.0 (all-transition panel)
        result = compute_titv(variants)
        assert result == 999.0

    def test_all_transversions(self):
        variants = self._make_variants([("A", "C"), ("A", "T"), ("G", "C"), ("G", "T")])
        # 0 Ti, 4 Tv → ratio = 0
        result = compute_titv(variants)
        assert result == 0.0

    def test_mixed(self):
        variants = self._make_variants([
            ("A", "G"), ("C", "T"),    # 2 Ti
            ("A", "C"), ("G", "T"),    # 2 Tv
        ])
        assert compute_titv(variants) == pytest.approx(1.0, abs=0.01)

    def test_skip_indels(self):
        variants = [
            {"ref": "A", "alt": "G", "variant_type": "SNV"},
            {"ref": "ATG", "alt": "A", "variant_type": "INDEL"},  # should be skipped
        ]
        # 1 Ti only (A->G), no Tv → sentinel 999.0
        assert compute_titv(variants) == 999.0

    def test_empty_list(self):
        assert compute_titv([]) == 0.0

    def test_typical_wes_ratio(self):
        """A ~2.5 Ti/Tv is typical for WES SNVs."""
        ti = [("A", "G")] * 25
        tv = [("A", "C")] * 10
        variants = self._make_variants(ti + tv)
        ratio = compute_titv(variants)
        assert ratio == pytest.approx(2.5, abs=0.01)


# ---------------------------------------------------------------------------
# extract_vep_csq
# ---------------------------------------------------------------------------

class TestExtractVepCsq:
    FIELDS = [
        "Allele", "Gene", "SYMBOL", "Consequence",
        "HGVSc", "HGVSp", "gnomAD_AF", "CADD_phred", "ClinVar_CLNSIG",
    ]

    def _csq(self, **kwargs):
        defaults = {
            "Allele": "G", "Gene": "ENSG001", "SYMBOL": "BRCA1",
            "Consequence": "missense_variant", "HGVSc": "NM_007294.4:c.100A>G",
            "HGVSp": "NP_009225.1:p.Thr34Ala", "gnomAD_AF": "0.0001",
            "CADD_phred": "28.5", "ClinVar_CLNSIG": "Pathogenic",
        }
        defaults.update(kwargs)
        return "|".join(defaults[f] for f in self.FIELDS)

    def test_gene_symbol_extracted(self):
        result = extract_vep_csq(self._csq(), self.FIELDS)
        assert result["gene"] == "BRCA1"

    def test_consequence_extracted(self):
        result = extract_vep_csq(self._csq(), self.FIELDS)
        assert result["consequence"] == "missense_variant"

    def test_hgvs_c_stripped_of_transcript(self):
        result = extract_vep_csq(self._csq(), self.FIELDS)
        # Should strip 'NM_007294.4:' prefix
        assert ":" not in result["hgvs_c"] or result["hgvs_c"].startswith("c.")

    def test_gnomad_af_parsed(self):
        result = extract_vep_csq(self._csq(), self.FIELDS)
        assert result["gnomad_af"] == pytest.approx(0.0001, rel=1e-3)

    def test_cadd_phred_parsed(self):
        result = extract_vep_csq(self._csq(), self.FIELDS)
        assert result["cadd_phred"] == pytest.approx(28.5, abs=0.01)

    def test_clinvar_extracted(self):
        result = extract_vep_csq(self._csq(), self.FIELDS)
        assert result["clinvar"] == "Pathogenic"

    def test_missing_gene_returns_empty(self):
        csq = "|".join(["G", "", "", "intergenic_variant", "", "", "0", "5", ""])
        result = extract_vep_csq(csq, self.FIELDS)
        assert result["gene"] == ""

    def test_multiple_transcripts_picks_most_severe(self):
        """With two transcripts, the more severe consequence should win."""
        tx1 = self._csq(Consequence="synonymous_variant", CADD_phred="5")
        tx2 = self._csq(SYMBOL="TP53", Consequence="stop_gained", CADD_phred="38")
        csq_str = f"{tx1},{tx2}"
        result = extract_vep_csq(csq_str, self.FIELDS)
        assert result["consequence"] == "stop_gained"

    def test_empty_af_returns_zero(self):
        result = extract_vep_csq(self._csq(gnomAD_AF="."), self.FIELDS)
        assert result["gnomad_af"] == 0.0

    def test_hgvsp_percent_decoded(self):
        result = extract_vep_csq(
            self._csq(HGVSp="NP_009225.1:p.Leu100%3D"),
            self.FIELDS,
        )
        assert "=" in result["hgvs_p"]


# ---------------------------------------------------------------------------
# _parse_vep_format_string
# ---------------------------------------------------------------------------

class TestParseVepFormatString:
    def test_parses_format_key(self):
        desc = 'Consequence annotations. Format: Allele|Gene|SYMBOL|Consequence'
        fields = _parse_vep_format_string(desc)
        assert fields == ["Allele", "Gene", "SYMBOL", "Consequence"]

    def test_returns_empty_when_no_format(self):
        desc = "Some description without a Format key"
        assert _parse_vep_format_string(desc) == []

    def test_strips_quotes(self):
        desc = 'Consequence annotations. Format: "Allele|Gene"'
        fields = _parse_vep_format_string(desc)
        assert "Allele" in fields


# ---------------------------------------------------------------------------
# _is_more_severe
# ---------------------------------------------------------------------------

class TestIsMorSevere:
    def test_frameshift_more_severe_than_missense(self):
        assert _is_more_severe("frameshift_variant", "missense_variant")

    def test_stop_gained_more_severe_than_synonymous(self):
        assert _is_more_severe("stop_gained", "synonymous_variant")

    def test_synonymous_not_more_severe_than_missense(self):
        assert not _is_more_severe("synonymous_variant", "missense_variant")

    def test_unknown_consequence_not_more_severe(self):
        assert not _is_more_severe("unknown_effect", "missense_variant")
