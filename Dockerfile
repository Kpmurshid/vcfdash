# ─────────────────────────────────────────────────────────────────────────────
# vcfdash — self-contained Docker image
#
# Bundles:
#   • vcfdash     (this tool)
#   • mosdepth    0.3.8   (coverage computation)
#   • Python      3.11    (slim base)
#   • htslib/pysam        (BAM/VCF I/O)
#
# Build:
#   docker build -t kpmurshid/vcfdash:latest .
#
# Run:
#   docker run --rm -v $(pwd):/data kpmurshid/vcfdash:latest \
#     python -m vcfdash \
#       --bam  /data/sample.bam \
#       --vcf  /data/sample.vep.vcf.gz \
#       --bed  /data/panel.bed \
#       --sample-id SAMPLE01 \
#       --outdir /data/reports
#
# Pull as Singularity SIF (HPC):
#   singularity pull vcfdash.sif docker://kpmurshid/vcfdash:latest
#   singularity exec --bind /data vcfdash.sif python -m vcfdash ...
# ─────────────────────────────────────────────────────────────────────────────

FROM python:3.11-slim

LABEL maintainer="Kpmurshid <https://github.com/Kpmurshid>"
LABEL org.opencontainers.image.title="vcfdash"
LABEL org.opencontainers.image.description="Clinical variant and coverage dashboard for WES/WGS"
LABEL org.opencontainers.image.source="https://github.com/Kpmurshid/vcfdash"
LABEL org.opencontainers.image.licenses="MIT"

# ── System dependencies ───────────────────────────────────────────────────────
# curl/xz-utils: download mosdepth static binary
# libz-dev, libbz2-dev, liblzma-dev, libcurl4: required by pysam/htslib
RUN apt-get update && apt-get install -y --no-install-recommends \
        curl \
        xz-utils \
        libz-dev \
        libbz2-dev \
        liblzma-dev \
        libcurl4 \
        libssl-dev \
        procps \
    && rm -rf /var/lib/apt/lists/*

# ── Install mosdepth static binary ───────────────────────────────────────────
# mosdepth releases a fully static x86_64 binary — no dependencies needed.
ARG MOSDEPTH_VERSION=0.3.8
RUN curl -fsSL \
    "https://github.com/brentp/mosdepth/releases/download/v${MOSDEPTH_VERSION}/mosdepth" \
    -o /usr/local/bin/mosdepth \
    && chmod +x /usr/local/bin/mosdepth \
    && mosdepth --version

# ── Install vcfdash Python package ───────────────────────────────────────────
WORKDIR /opt/vcfdash

# Copy all files needed by hatchling to build the package
COPY pyproject.toml README.md ./
COPY vcfdash/ ./vcfdash/

# Install Python dependencies and vcfdash itself (non-editable for container use)
RUN pip install --no-cache-dir --upgrade pip \
    && pip install --no-cache-dir \
        jinja2>=3.0 \
        pysam>=0.22 \
    && pip install --no-cache-dir .

# ── Pre-fetch D3.js so the image works fully offline ─────────────────────────
RUN python -c "from vcfdash.report import _get_d3_js; _get_d3_js()" \
    && echo "D3.js cached."

# ── Create working directory for user data ───────────────────────────────────
RUN mkdir -p /data
WORKDIR /data

# ── Sanity check ─────────────────────────────────────────────────────────────
RUN python -m vcfdash --version \
    && mosdepth --version

# ── Default entry point ───────────────────────────────────────────────────────
# Pass --help to see all options, or override the full command:
#   docker run kpmurshid/vcfdash python -m vcfdash --bam ...
ENTRYPOINT ["python", "-m", "vcfdash"]
CMD ["--help"]
