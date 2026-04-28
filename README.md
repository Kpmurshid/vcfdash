# vcfdash

**Lightweight, pipeline-native clinical variant and coverage dashboard.**

vcfdash takes a sorted BAM, an annotated VCF, and a target BED file for a single WES/WGS sample and produces a **fully self-contained HTML report** — one file, no server, no internet required — plus an optional machine-readable JSON summary.

[![Python ≥3.10](https://img.shields.io/badge/python-≥3.10-blue)](https://www.python.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Tests](https://img.shields.io/badge/tests-87%20passed-brightgreen)]()

> **GitHub:** [https://github.com/Kpmurshid/vcfdash](https://github.com/Kpmurshid/vcfdash)

---

## Features

| Feature | Details |
|---------|---------|
| **Zero-dependency output** | Single `.html` file opens in any modern browser — all CSS, JS, and data inlined |
| **Three interactive panels** | ① QC summary · ② Coverage table · ③ Variant table |
| **VEP & SnpEff support** | Auto-detects `CSQ` or `ANN` annotation fields; graceful degradation if absent |
| **D3.js sparklines** | Per-variant 200 bp coverage bar chart rendered on row expand |
| **Sortable & filterable tables** | Gene/region search, status filter, ClinVar filter, free-text search |
| **JSON export** | Machine-readable summary for downstream programmatic use |
| **Nextflow module** | nf-core-compatible process module included |
| **Singularity / Apptainer support** | Run mosdepth from a SIF image — no root, no conda needed on HPC |
| **Light-themed, full-width UI** | Responsive layout that fills the screen at any resolution |

---

## Installation

### pip (recommended)

```bash
pip install vcfdash
```

### From source

```bash
git clone https://github.com/Kpmurshid/vcfdash.git
cd vcfdash
pip install -e ".[dev]"
```

### Dependencies

| Tool | Version | Purpose |
|------|---------|---------|
| Python | ≥ 3.10 | Runtime |
| [mosdepth](https://github.com/brentp/mosdepth) | ≥ 0.3.3 | Coverage computation |
| jinja2 | ≥ 3.0 | HTML templating |
| pysam | ≥ 0.22 | BAM/VCF reading |

---

## mosdepth — Installation Options

vcfdash calls mosdepth internally. You can provide it in three ways:

### Option 1 — conda / mamba (recommended for workstations)

```bash
conda install -c bioconda mosdepth
```

### Option 2 — Singularity / Apptainer SIF (recommended for HPC)

If mosdepth is **not in PATH**, vcfdash auto-detects SIF images in common locations.
You can also supply the SIF path explicitly:

```bash
# Build or pull a SIF image
singularity pull mosdepth.sif docker://quay.io/biocontainers/mosdepth:0.3.8--hd299d5a_0

# Then pass it to vcfdash
python -m vcfdash \
  --bam  sample.bam \
  --vcf  sample.vep.vcf.gz \
  --bed  panel.bed \
  --sample-id SAMPLE01 \
  --mosdepth-sif /path/to/mosdepth.sif
```

vcfdash will run:
```
singularity exec --bind <dirs> /path/to/mosdepth.sif mosdepth <args>
```

Auto-detected SIF locations (checked in order):
- `/COLD_STORAGE/software/tools/mosdepth/mosdepth.sif`
- `/opt/software/mosdepth/mosdepth.sif`
- `~/singularity/mosdepth.sif`

### Option 3 — Docker

```bash
docker pull quay.io/biocontainers/mosdepth:0.3.8--hd299d5a_0
```

Use a wrapper script that mounts directories and delegates to Docker.

---

## Quick Start

```bash
python -m vcfdash \
  --bam  sample.bam \
  --vcf  sample.vep.vcf.gz \
  --bed  panel.bed \
  --sample-id SAMPLE01 \
  --outdir reports/ \
  --genome hg38 \
  --json \
  --no-sparklines
```

Output:
```
reports/SAMPLE01_report.html   ← Self-contained HTML (open in browser)
reports/SAMPLE01_report.json   ← Machine-readable JSON summary
```

---

## Full Usage

```
usage: vcfdash [-h] [--version]
               --bam PATH --vcf PATH --bed PATH --sample-id STRING
               [--outdir PATH] [--genome STRING]
               [--min-dp INT] [--min-dp2 INT]
               [--config PATH] [--json] [--no-sparklines]
               [--mosdepth-sif PATH]

required arguments:
  --bam PATH          Sorted, indexed BAM (.bam + .bai must exist)
  --vcf PATH          Annotated VCF (.vcf or .vcf.gz; VEP/SnpEff preferred)
  --bed PATH          Target regions BED (exome capture or gene panel)
  --sample-id STRING  Sample identifier used in report title and filenames

optional arguments:
  --outdir PATH       Output directory (default: current directory)
  --genome STRING     hg38 (default) | hg19 | GRCh38 | GRCh37
  --min-dp INT        Primary depth threshold for PASS/WARN/FAIL (default: 20)
  --min-dp2 INT       Secondary depth threshold (default: 10)
  --config PATH       JSON config file for custom thresholds
  --json              Also write a machine-readable JSON summary
  --no-sparklines     Skip per-variant sparklines — faster, smaller HTML
  --mosdepth-sif PATH Path to mosdepth Singularity/Apptainer SIF image
```

---

## Output

### ① QC Summary panel

Key metrics displayed as colour-coded cards (green / amber / red):

| Metric | Description |
|--------|-------------|
| **Mean Coverage (on-BED)** | Length-weighted mean depth over supplied BED intervals |
| **Bases ≥20x** | % of BED bases at or above primary depth threshold |
| **Bases ≥10x** | % of BED bases at or above secondary depth threshold |
| **Total Variants** | Count of all variants (PASS + flagged) |
| **SNV / INDEL / MNV** | Variant type breakdown |
| **Ti/Tv ratio** | Transition/Transversion ratio (expected 2.0–3.3 for WES) |

> **Note on coverage:** The 108x shown for a Twist WES sample is **mean bait/capture coverage** (over the 37.5 Mb Twist BED). Picard HsMetrics reports a lower `MEAN_TARGET_COVERAGE` (e.g. 64.5x) because it uses a smaller, merged exon-only BED. Both are correct — they measure different territories.

### ② Coverage panel

Sortable, searchable table. Each row = one chromosome (for exome-capture BEDs) or one gene/region (for gene-panel BEDs with named intervals). Click a row to jump to that gene's variants.

| Column | Description |
|--------|-------------|
| Region / Gene | Chromosome or gene name from BED column 4 |
| Mean Depth | Length-weighted mean coverage |
| ≥20x % | Estimated fraction of bases above primary threshold |
| ≥10x % | Estimated fraction above secondary threshold |
| Status | ✅ PASS · ⚠️ WARN · ❌ FAIL |

### ③ Variant panel

Filterable by gene, consequence, ClinVar classification, or free text. Loads in batches of 200. Each row expands to show:
- Full VEP annotation (HGVSc, HGVSp, consequence, transcript)
- Population frequency (gnomAD AF)
- ClinVar classification
- CADD phred score
- Genotype fields (GT, DP, GQ, AD, VAF)
- 200 bp coverage sparkline (if `--no-sparklines` not set)

### JSON export (`--json`)

```json
{
  "sample_id": "SAMPLE01",
  "genome": "hg38",
  "generated_at": "2025-01-01T00:00:00Z",
  "thresholds": { "min_dp": 20, "min_dp2": 10 },
  "qc_summary": { "mean_coverage": 108.0, "pct_bases_above_20x": 95.8, ... },
  "gene_coverage": [ { "gene": "chr1", "mean_depth": 108.8, "status": "PASS", ... } ],
  "variants": [ { "chrom": "chr1", "pos": 925952, "ref": "G", "alt": "A", ... } ]
}
```

---

## Nextflow Integration

Copy the module into your pipeline:

```bash
cp modules/nextflow/vcfdash.nf <your-pipeline>/modules/vcfdash.nf
```

```nextflow
include { VCFDASH } from './modules/vcfdash.nf'

workflow {
    VCFDASH(
        tuple(meta, bam, bai, vcf, tbi, bed)
    )
}
```

Pipeline parameters:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `params.vcfdash_json`    | `false` | Enable JSON export |
| `params.genome`          | `hg38`  | Reference genome build |
| `params.vcfdash_min_dp`  | `20`    | Primary depth threshold |
| `params.vcfdash_min_dp2` | `10`    | Secondary depth threshold |
| `params.mosdepth_sif`    | `null`  | Path to mosdepth SIF image (HPC) |

---

## Configuration file

Create a JSON config to override defaults:

```json
{
  "pass_pct_20x": 95.0,
  "warn_pct_20x": 80.0,
  "pass_depth_fraction": 1.0,
  "warn_depth_fraction": 0.5
}
```

Pass with `--config my_config.json`. CLI flags always take precedence over the config file.

---

## VCF Annotation Compatibility

| Tool | INFO field | Support |
|------|-----------|---------|
| Ensembl VEP | `CSQ` | ✅ Primary |
| SnpEff | `ANN` | ✅ Fallback |
| ANNOVAR | INFO fields | Partial (best via VEP post-processing) |
| Unannotated VCF | — | ✅ Raw variant data shown; warning in report header |

---

## Running Tests

```bash
pip install -e ".[dev]"
pytest tests/ -v
# 87 passed in < 1 s
```

---

## Error Handling

| Situation | Behaviour |
|-----------|-----------|
| BAM index missing | `ERROR: BAM index not found. Run: samtools index {bam}` |
| mosdepth not in PATH and no SIF | `ERROR: mosdepth not found.` with install hint |
| SIF not executable | `ERROR: SIF image not found or not readable: {path}` |
| No CSQ/ANN field | Warning banner in report header; raw variant positions shown |
| BED without column 4 name | `chrom:start-end` used as region label |
| No AD field in VCF | VAF shown as N/A; sparkline skipped |

All errors write to stderr. vcfdash never silently fails.

---

## Project Structure

```
vcfdash/
├── vcfdash/
│   ├── __init__.py        # Version
│   ├── __main__.py        # python -m vcfdash entry point
│   ├── cli.py             # Argument parsing + pipeline orchestration
│   ├── coverage.py        # mosdepth runner + BED/summary parser
│   ├── variants.py        # VCF parser + VEP/SnpEff annotation extractor
│   ├── report.py          # HTML renderer + JSON exporter
│   ├── utils.py           # Shared helpers (file checks, config loader)
│   └── templates/
│       ├── report.html.j2        # Jinja2 HTML template
│       └── assets/
│           ├── style.css         # Light-theme CSS (inlined into HTML)
│           ├── viz.js            # D3 sparklines + table sort/filter
│           └── d3.min.js         # D3.js v7 (auto-downloaded on first run)
├── modules/
│   └── nextflow/vcfdash.nf       # nf-core-compatible Nextflow module
├── tests/
│   ├── conftest.py
│   ├── test_coverage.py
│   ├── test_variants.py
│   └── test_report.py
├── examples/
│   └── smoke_test.py
├── pyproject.toml
└── README.md
```

---

## Contributing

Contributions are welcome! Please open an issue or pull request at [github.com/Kpmurshid/vcfdash](https://github.com/Kpmurshid/vcfdash).

1. Fork the repository
2. Install dev dependencies: `pip install -e ".[dev]"`
3. Run tests: `pytest`
4. Submit a PR against `main`

---

## License

MIT — see [LICENSE](LICENSE).

---

## Acknowledgements

- [mosdepth](https://github.com/brentp/mosdepth) — fast coverage computation
- [D3.js](https://d3js.org/) — sparkline visualisations
- [Jinja2](https://jinja.palletsprojects.com/) — HTML templating
- [Ensembl VEP](https://www.ensembl.org/vep) — variant annotation
