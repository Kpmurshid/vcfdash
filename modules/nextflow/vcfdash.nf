/*
 * vcfdash Nextflow process module
 * nf-core compatible
 *
 * Produces a self-contained HTML report (and optional JSON) from
 * BAM + VCF + BED inputs for a single sample.
 *
 * Required parameters (via task.ext.args or params):
 *   params.vcfdash_json  (bool)  — enable JSON export (default: false)
 *
 * Meta map keys used:
 *   meta.id  — sample identifier passed to --sample-id
 */

process VCFDASH {
    tag "${meta.id}"
    label 'process_low'

    conda 'bioconda::vcfdash=0.1.0'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vcfdash:0.1.0--pyhdfd78af_0' :
        'quay.io/biocontainers/vcfdash:0.1.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai), path(vcf), path(tbi), path(bed)

    output:
    tuple val(meta), path("*.html"),            emit: report
    tuple val(meta), path("*.json"),            emit: json,    optional: true
    path  "versions.yml",                       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args   ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def json_flag   = params.vcfdash_json ? '--json' : ''
    def genome      = params.genome       ?: 'hg38'
    def min_dp      = params.vcfdash_min_dp  ?: 20
    def min_dp2     = params.vcfdash_min_dp2 ?: 10
    """
    vcfdash \\
        --bam ${bam} \\
        --vcf ${vcf} \\
        --bed ${bed} \\
        --sample-id ${prefix} \\
        --outdir . \\
        --genome ${genome} \\
        --min-dp ${min_dp} \\
        --min-dp2 ${min_dp2} \\
        ${json_flag} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcfdash: \$(vcfdash --version 2>&1 | sed 's/vcfdash //')
        mosdepth: \$(mosdepth --version 2>&1 | sed 's/mosdepth //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_report.html
    touch ${prefix}_report.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcfdash: \$(vcfdash --version 2>&1 | sed 's/vcfdash //')
        mosdepth: \$(mosdepth --version 2>&1 | sed 's/mosdepth //')
    END_VERSIONS
    """
}
