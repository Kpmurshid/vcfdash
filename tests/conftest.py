"""
conftest.py — Shared pytest fixtures for vcfdash tests.
"""

from __future__ import annotations

import gzip
import io
import json
import textwrap
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import pytest


# ---------------------------------------------------------------------------
# Minimal synthetic args namespace
# ---------------------------------------------------------------------------

@pytest.fixture
def args(tmp_path):
    """Return a minimal argparse-like namespace for testing."""
    return SimpleNamespace(
        bam       = str(tmp_path / "sample.bam"),
        vcf       = str(tmp_path / "sample.vcf"),
        bed       = str(tmp_path / "panel.bed"),
        sample_id = "TEST_SAMPLE",
        outdir    = str(tmp_path / "out"),
        genome    = "hg38",
        min_dp    = 20,
        min_dp2   = 10,
        json      = False,
        no_sparklines = True,
        config    = None,
    )


# ---------------------------------------------------------------------------
# Synthetic BED file
# ---------------------------------------------------------------------------

@pytest.fixture
def bed_file(tmp_path) -> Path:
    """Write a minimal 4-column BED file."""
    path = tmp_path / "panel.bed"
    path.write_text(
        "chr1\t100\t200\tGENE_A\n"
        "chr1\t300\t400\tGENE_B\n"
        "chr2\t500\t600\tGENE_C\n",
        encoding="utf-8",
    )
    return path


# ---------------------------------------------------------------------------
# Synthetic mosdepth output files
# ---------------------------------------------------------------------------

@pytest.fixture
def mosdepth_prefix(tmp_path) -> str:
    """
    Write synthetic mosdepth output files under tmp_path/mosdepth/sample
    and return the prefix string.
    """
    md_dir = tmp_path / "mosdepth"
    md_dir.mkdir()
    prefix = str(md_dir / "sample")

    # regions.bed.gz
    regions_content = (
        "chr1\t100\t200\tGENE_A\t85.3\n"
        "chr1\t300\t400\tGENE_B\t12.1\n"
        "chr2\t500\t600\tGENE_C\t7.4\n"
    )
    with gzip.open(f"{prefix}.regions.bed.gz", "wt") as fh:
        fh.write(regions_content)

    # mosdepth.summary.txt
    summary_content = (
        "chrom\tlength\tbases\tmean\tmin\tmax\n"
        "chr1\t200\t17060\t85.30\t30\t120\n"
        "chr2\t100\t740\t7.40\t0\t15\n"
        "total\t300\t17800\t59.33\t0\t120\n"
    )
    Path(f"{prefix}.mosdepth.summary.txt").write_text(summary_content)

    return prefix


@pytest.fixture
def coverage_data(mosdepth_prefix) -> dict[str, Any]:
    """Return parsed coverage data from synthetic mosdepth output."""
    from vcfdash.coverage import parse_mosdepth_output
    return parse_mosdepth_output(
        prefix       = mosdepth_prefix,
        min_dp       = 20,
        min_dp2      = 10,
        no_sparklines= True,
    )


# ---------------------------------------------------------------------------
# Synthetic VCF content helpers
# ---------------------------------------------------------------------------

VCF_HEADER_VEP = textwrap.dedent("""\
    ##fileformat=VCFv4.2
    ##FILTER=<ID=PASS,Description="All filters passed">
    ##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Gene|SYMBOL|Consequence|HGVSc|HGVSp|gnomAD_AF|CADD_phred|ClinVar_CLNSIG">
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
    ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
    ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
""")

VCF_HEADER_PLAIN = textwrap.dedent("""\
    ##fileformat=VCFv4.2
    ##FILTER=<ID=PASS,Description="All filters passed">
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
    ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
    ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
""")


@pytest.fixture
def vcf_vep(tmp_path) -> Path:
    """Minimal VEP-annotated VCF with two variants."""
    path = tmp_path / "sample.vep.vcf"
    # CSQ format: Allele|Gene|SYMBOL|Consequence|HGVSc|HGVSp|gnomAD_AF|CADD_phred|ClinVar_CLNSIG
    content = VCF_HEADER_VEP + (
        "chr1\t150\t.\tA\tG\t.\tPASS\t"
        "CSQ=G|ENSG001|GENE_A|missense_variant|c.100A>G|p.Thr34Ala|0.0001|28.5|Uncertain_significance\t"
        "GT:DP:GQ:AD\t0/1:45:99:22,23\n"
        "chr1\t350\t.\tATG\tA\t.\tLowQual\t"
        "CSQ=A|ENSG002|GENE_B|frameshift_variant|c.250_252del|p.Leu84del|.|35.2|\t"
        "GT:DP:GQ:AD\t0/1:18:45:10,8\n"
    )
    path.write_text(content, encoding="utf-8")
    return path


@pytest.fixture
def vcf_plain(tmp_path) -> Path:
    """Plain VCF without annotation fields."""
    path = tmp_path / "sample.plain.vcf"
    content = VCF_HEADER_PLAIN + (
        "chr2\t550\t.\tC\tT\t.\tPASS\t.\tGT:DP:GQ:AD\t1/1:60:99:0,60\n"
    )
    path.write_text(content, encoding="utf-8")
    return path


# ---------------------------------------------------------------------------
# Synthetic variant data
# ---------------------------------------------------------------------------

@pytest.fixture
def variant_data() -> list[dict]:
    """Return a small synthetic list of variant dicts."""
    return [
        {
            "chrom": "chr1", "pos": 150, "ref": "A", "alt": "G",
            "filter": "PASS", "gt": "0/1", "dp": 45, "gq": 99,
            "vaf": 0.511, "ad": [22, 23], "variant_type": "SNV",
            "gene": "GENE_A", "consequence": "missense_variant",
            "hgvs_c": "c.100A>G", "hgvs_p": "p.Thr34Ala",
            "gnomad_af": 0.0001, "cadd_phred": 28.5,
            "clinvar": "Uncertain_significance",
        },
        {
            "chrom": "chr1", "pos": 350, "ref": "ATG", "alt": "A",
            "filter": "LowQual", "gt": "0/1", "dp": 18, "gq": 45,
            "vaf": 0.444, "ad": [10, 8], "variant_type": "INDEL",
            "gene": "GENE_B", "consequence": "frameshift_variant",
            "hgvs_c": "c.250_252del", "hgvs_p": "p.Leu84del",
            "gnomad_af": 0.0, "cadd_phred": 35.2, "clinvar": "",
        },
        {
            "chrom": "chr2", "pos": 550, "ref": "C", "alt": "T",
            "filter": "PASS", "gt": "1/1", "dp": 60, "gq": 99,
            "vaf": 1.0, "ad": [0, 60], "variant_type": "SNV",
            "gene": "GENE_C", "consequence": "synonymous_variant",
            "hgvs_c": "c.300C>T", "hgvs_p": "p.Leu100=",
            "gnomad_af": 0.05, "cadd_phred": 5.1, "clinvar": "Benign",
        },
    ]
