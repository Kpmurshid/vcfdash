"""
variants.py — VCF parsing, INFO field extraction, VAF calculation.

Primary parser: cyvcf2
Fallback:       pysam VariantFile
Annotation:     VEP CSQ field (preferred) → SnpEff ANN field → raw data only
"""

from __future__ import annotations

import sys
import warnings
from typing import Any

from .utils import safe_float, safe_int

# ---------------------------------------------------------------------------
# VEP CSQ field names we care about (in order of preference)
# ---------------------------------------------------------------------------

VEP_FIELDS_OF_INTEREST = [
    "Gene",
    "SYMBOL",
    "Consequence",
    "HGVSc",
    "HGVSp",
    "gnomAD_AF",
    "gnomADe_AF",
    "gnomADg_AF",
    "CADD_phred",
    "CADD_PHRED",
    "ClinVar_CLNSIG",
    "CLIN_SIG",
    "CLINVAR_CLNSIG",
    "Existing_variation",
    "Feature",
]

# SnpEff ANN sub-fields (pipe-delimited)
SNPEFF_ANN_FIELDS = [
    "Allele", "Annotation", "Annotation_Impact", "Gene_Name", "Gene_ID",
    "Feature_Type", "Feature_ID", "Transcript_BioType", "Rank", "HGVS.c",
    "HGVS.p", "cDNA.pos / cDNA.length", "CDS.pos / CDS.length",
    "AA.pos / AA.length", "Distance", "ERRORS / WARNINGS / INFO",
]


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def parse_vcf(vcf_path: str, no_sparklines: bool = False) -> tuple[list[dict], list[str]]:
    """
    Parse *vcf_path* and return (variant_records, warnings).

    Each variant record is a dict with keys:
      chrom, pos, ref, alt, filter, gt, dp, gq, vaf,
      gene, consequence, hgvs_c, hgvs_p, gnomad_af,
      cadd_phred, clinvar, variant_type, ad
    """
    warnings_list: list[str] = []

    try:
        from cyvcf2 import VCF  # type: ignore
        return _parse_with_cyvcf2(vcf_path, warnings_list)
    except ImportError:
        pass

    try:
        import pysam  # type: ignore
        return _parse_with_pysam(vcf_path, warnings_list)
    except ImportError:
        pass

    from .utils import _fatal
    _fatal(
        "Neither cyvcf2 nor pysam is installed.\n"
        "Install via:  pip install cyvcf2\n"
        "or:           pip install pysam"
    )
    return [], []  # unreachable


# ---------------------------------------------------------------------------
# cyvcf2 implementation
# ---------------------------------------------------------------------------

def _parse_with_cyvcf2(vcf_path: str, warnings_list: list[str]) -> tuple[list[dict], list[str]]:
    from cyvcf2 import VCF  # type: ignore

    vcf = VCF(vcf_path)

    # Detect annotation style
    csq_fields  = _detect_vep_csq_fields(vcf)
    ann_fields  = _detect_snpeff_ann_fields(vcf)
    annotation_style: str

    if csq_fields:
        annotation_style = "vep"
    elif ann_fields:
        annotation_style = "snpeff"
    else:
        annotation_style = "none"
        warnings_list.append(
            "No CSQ (VEP) or ANN (SnpEff) annotation field found in VCF. "
            "Displaying raw variant data only."
        )

    records: list[dict] = []

    for variant in vcf:
        rec = _extract_cyvcf2_record(variant, csq_fields, ann_fields, annotation_style)
        records.append(rec)

    vcf.close()
    return records, warnings_list


def _extract_cyvcf2_record(
    variant,
    csq_fields: list[str],
    ann_fields: list[str],
    annotation_style: str,
) -> dict[str, Any]:
    """Extract a single variant record from a cyvcf2 variant object."""
    chrom  = variant.CHROM
    pos    = variant.POS
    ref    = variant.REF
    alts   = variant.ALT
    alt    = alts[0] if alts else "."
    filt   = _format_filter(variant.FILTER)

    vtype  = classify_variant_type(ref, alt)

    # FORMAT fields — GT, DP, GQ, AD
    gt_str = "."
    dp     = 0
    gq     = 0
    vaf    = None
    ad     = None

    try:
        gts = variant.genotypes  # list of [allele1, allele2, phased]
        if gts:
            g = gts[0]
            a1, a2 = g[0], g[1]
            phased = g[2]
            sep    = "|" if phased else "/"
            gt_str = f"{a1}{sep}{a2}"
    except Exception:
        pass

    try:
        dp_val = variant.format("DP")
        if dp_val is not None and len(dp_val) > 0:
            dp = safe_int(dp_val[0][0])
    except Exception:
        pass

    try:
        gq_val = variant.format("GQ")
        if gq_val is not None and len(gq_val) > 0:
            gq = safe_int(gq_val[0][0])
    except Exception:
        pass

    try:
        ad_val = variant.format("AD")
        if ad_val is not None and len(ad_val) > 0:
            ad_list = list(ad_val[0])
            ad = ad_list
            ref_ad = safe_int(ad_list[0])
            alt_ad = safe_int(ad_list[1]) if len(ad_list) > 1 else 0
            total  = ref_ad + alt_ad
            vaf    = round(alt_ad / total, 4) if total > 0 else None
    except Exception:
        pass

    # Annotation extraction
    anno = _empty_annotation()

    if annotation_style == "vep":
        csq_raw = variant.INFO.get("CSQ")
        if csq_raw:
            anno = extract_vep_csq(csq_raw, csq_fields)
    elif annotation_style == "snpeff":
        ann_raw = variant.INFO.get("ANN")
        if ann_raw:
            anno = _extract_snpeff_ann(ann_raw, ann_fields)

    return {
        "chrom":       chrom,
        "pos":         pos,
        "ref":         ref,
        "alt":         alt,
        "filter":      filt,
        "gt":          gt_str,
        "dp":          dp,
        "gq":          gq,
        "vaf":         vaf,
        "ad":          ad,
        "variant_type": vtype,
        **anno,
    }


# ---------------------------------------------------------------------------
# pysam fallback implementation
# ---------------------------------------------------------------------------

def _parse_with_pysam(vcf_path: str, warnings_list: list[str]) -> tuple[list[dict], list[str]]:
    import pysam  # type: ignore

    vcf = pysam.VariantFile(vcf_path)

    # Detect annotation style from header
    header_info = dict(vcf.header.info)
    csq_fields: list[str] = []
    ann_fields: list[str] = []
    annotation_style = "none"

    if "CSQ" in header_info:
        desc = header_info["CSQ"].description or ""
        csq_fields = _parse_vep_format_string(desc)
        annotation_style = "vep"
    elif "ANN" in header_info:
        annotation_style = "snpeff"
        ann_fields = SNPEFF_ANN_FIELDS
    else:
        warnings_list.append(
            "No CSQ (VEP) or ANN (SnpEff) annotation field found in VCF. "
            "Displaying raw variant data only."
        )

    records: list[dict] = []

    for variant in vcf.fetch():
        rec = _extract_pysam_record(variant, csq_fields, ann_fields, annotation_style)
        records.append(rec)

    vcf.close()
    return records, warnings_list


def _extract_pysam_record(
    variant,
    csq_fields: list[str],
    ann_fields: list[str],
    annotation_style: str,
) -> dict[str, Any]:
    """Extract a single variant record from a pysam variant."""
    chrom = variant.chrom
    pos   = variant.pos + 1   # pysam is 0-based
    ref   = variant.ref
    alts  = variant.alts or (".",)
    alt   = alts[0]
    filt  = ",".join(list(variant.filter)) if variant.filter else "PASS"
    vtype = classify_variant_type(ref, alt)

    gt_str = "."
    dp     = 0
    gq     = 0
    vaf    = None
    ad     = None

    try:
        sample = list(variant.samples.values())[0]
        gt_alleles = sample["GT"]
        if gt_alleles:
            gt_str = "/".join(str(a) if a is not None else "." for a in gt_alleles)
        dp  = safe_int(sample.get("DP", 0))
        gq  = safe_int(sample.get("GQ", 0))
        ad_val = sample.get("AD")
        if ad_val is not None:
            ad = list(ad_val)
            ref_ad = safe_int(ad[0])
            alt_ad = safe_int(ad[1]) if len(ad) > 1 else 0
            total  = ref_ad + alt_ad
            vaf    = round(alt_ad / total, 4) if total > 0 else None
    except Exception:
        pass

    anno = _empty_annotation()

    if annotation_style == "vep":
        csq_raw = variant.info.get("CSQ")
        if csq_raw:
            csq_str = csq_raw if isinstance(csq_raw, str) else ",".join(csq_raw)
            anno = extract_vep_csq(csq_str, csq_fields)
    elif annotation_style == "snpeff":
        ann_raw = variant.info.get("ANN")
        if ann_raw:
            ann_str = ann_raw if isinstance(ann_raw, str) else ",".join(ann_raw)
            anno = _extract_snpeff_ann(ann_str, ann_fields)

    return {
        "chrom":       chrom,
        "pos":         pos,
        "ref":         ref,
        "alt":         alt,
        "filter":      filt,
        "gt":          gt_str,
        "dp":          dp,
        "gq":          gq,
        "vaf":         vaf,
        "ad":          ad,
        "variant_type": vtype,
        **anno,
    }


# ---------------------------------------------------------------------------
# VEP CSQ extraction
# ---------------------------------------------------------------------------

def _detect_vep_csq_fields(vcf) -> list[str]:
    """Extract VEP CSQ field names from cyvcf2 VCF header."""
    try:
        desc = vcf.get_header_type("CSQ").get("Description", "")
        return _parse_vep_format_string(desc)
    except Exception:
        return []


def _parse_vep_format_string(description: str) -> list[str]:
    """
    Parse 'Format: Allele|Gene|Feature|...' from VEP CSQ Description.
    """
    if "Format:" in description:
        fmt_part = description.split("Format:")[-1].strip().strip('"').strip("'")
        return [f.strip() for f in fmt_part.split("|")]
    return []


def extract_vep_csq(csq_string: str, csq_fields: list[str]) -> dict[str, Any]:
    """
    Parse a VEP CSQ INFO value into structured annotation dict.

    CSQ is a comma-delimited list of transcripts, each pipe-delimited.
    We pick the first transcript with a canonical/protein-coding consequence.
    """
    best    = _empty_annotation()
    transcripts = csq_string.split(",")

    for tx_str in transcripts:
        tx_vals = tx_str.split("|")
        tx = dict(zip(csq_fields, tx_vals))

        # Prefer transcripts that have a gene symbol
        gene = tx.get("SYMBOL") or tx.get("Gene") or ""
        if not gene:
            continue

        consequence  = tx.get("Consequence", "")
        hgvs_c       = tx.get("HGVSc", "").split(":")[-1]
        hgvs_p       = tx.get("HGVSp", "").split(":")[-1].replace("%3D", "=")

        gnomad_af = 0.0
        for af_key in ("gnomAD_AF", "gnomADe_AF", "gnomADg_AF", "MAX_AF"):
            raw = tx.get(af_key, "")
            if raw and raw not in (".", ""):
                gnomad_af = safe_float(raw)
                break

        cadd_phred = safe_float(tx.get("CADD_phred") or tx.get("CADD_PHRED") or "0")

        clinvar = (
            tx.get("ClinVar_CLNSIG")
            or tx.get("CLIN_SIG")
            or tx.get("CLINVAR_CLNSIG")
            or ""
        )

        anno = {
            "gene":        gene,
            "consequence": consequence,
            "hgvs_c":      hgvs_c,
            "hgvs_p":      hgvs_p,
            "gnomad_af":   gnomad_af,
            "cadd_phred":  cadd_phred,
            "clinvar":     clinvar,
        }

        # Use this transcript if it's better than what we have
        if not best["gene"] or _is_more_severe(consequence, best["consequence"]):
            best = anno

    return best


# ---------------------------------------------------------------------------
# SnpEff ANN extraction
# ---------------------------------------------------------------------------

def _detect_snpeff_ann_fields(vcf) -> list[str]:
    """Check cyvcf2 header for ANN field."""
    try:
        vcf.get_header_type("ANN")
        return SNPEFF_ANN_FIELDS
    except Exception:
        return []


def _extract_snpeff_ann(ann_string: str, ann_fields: list[str]) -> dict[str, Any]:
    """
    Parse SnpEff ANN INFO value.

    ANN is a comma-delimited list; each element is pipe-delimited.
    """
    best = _empty_annotation()

    for ann_entry in ann_string.split(","):
        parts = ann_entry.split("|")
        ann   = dict(zip(ann_fields, parts))

        gene        = ann.get("Gene_Name", "")
        consequence = ann.get("Annotation", "")
        hgvs_c      = ann.get("HGVS.c", "")
        hgvs_p      = ann.get("HGVS.p", "")

        anno = {
            "gene":        gene,
            "consequence": consequence,
            "hgvs_c":      hgvs_c,
            "hgvs_p":      hgvs_p,
            "gnomad_af":   0.0,
            "cadd_phred":  0.0,
            "clinvar":     "",
        }

        if not best["gene"] or _is_more_severe(consequence, best["consequence"]):
            best = anno

    return best


# ---------------------------------------------------------------------------
# Variant classification and metrics
# ---------------------------------------------------------------------------

def classify_variant_type(ref: str, alt: str) -> str:
    """Return 'SNV', 'INDEL', or 'MNV'."""
    if alt in (".", "*"):
        return "SNV"
    if len(ref) == 1 and len(alt) == 1:
        return "SNV"
    if len(ref) != len(alt):
        return "INDEL"
    if len(ref) > 1 and len(ref) == len(alt):
        return "MNV"
    return "SNV"


def compute_titv(variants: list[dict]) -> float:
    """
    Compute Ti/Tv ratio from SNV variants.

    Transitions (Ti): A↔G, C↔T
    Transversions (Tv): A↔C, A↔T, G↔C, G↔T
    """
    ti_pairs = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}

    ti = 0
    tv = 0

    for v in variants:
        if v.get("variant_type") != "SNV":
            continue
        ref = (v.get("ref") or "").upper()
        alt = (v.get("alt") or "").upper()
        if len(ref) != 1 or len(alt) != 1:
            continue
        if (ref, alt) in ti_pairs:
            ti += 1
        else:
            tv += 1

    if ti == 0 and tv == 0:
        return 0.0
    if tv == 0:
        return 999.0   # all transitions — return a high sentinel value
    return round(ti / tv, 3)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _empty_annotation() -> dict[str, Any]:
    return {
        "gene":        "",
        "consequence": "",
        "hgvs_c":      "",
        "hgvs_p":      "",
        "gnomad_af":   0.0,
        "cadd_phred":  0.0,
        "clinvar":     "",
    }


def _format_filter(filt) -> str:
    """Format cyvcf2 FILTER field to string."""
    if filt is None:
        return "PASS"
    if isinstance(filt, (list, tuple)):
        return ",".join(str(f) for f in filt) if filt else "PASS"
    return str(filt)


# Severity order for consequence selection (higher index = more severe)
_SEVERITY_ORDER = [
    "intergenic_variant",
    "intron_variant",
    "synonymous_variant",
    "splice_region_variant",
    "missense_variant",
    "splice_donor_variant",
    "splice_acceptor_variant",
    "stop_gained",
    "frameshift_variant",
    "start_lost",
    "stop_lost",
    "transcript_ablation",
]


def _is_more_severe(consequence: str, current_best: str) -> bool:
    """Return True if *consequence* is more severe than *current_best*."""
    def rank(c: str) -> int:
        for i, term in enumerate(_SEVERITY_ORDER):
            if term in c:
                return i
        return -1

    return rank(consequence) > rank(current_best)
