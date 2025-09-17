process CUSTOM_BAMTOWIG {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda:pysam==0.18.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pysam:0.18.0--py39h5030a8b_2':
        'quay.io/biocontainers/pysam:0.18.0--py38h104f7d5_0' }"

    input:

    tuple val(meta), path(bam), path(bam_index), val(name), path(fasta)
    val transposon_type
    val ta_position
    
    output:
    tuple val(name), val(meta), path("*.wig"), emit: transit_wig
    path "versions.yml"           , emit: versions
    path "*.log"                  , emit: log

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    bam_to_wig.py \\
        $bam \\
        --outfile ${name}_${prefix}.wig \\
        --transposon $transposon_type \\
        --fasta $fasta \\
        --ta-position $ta_position \\
        $args \\
        > ${name}_${prefix}.bam_to_wig.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pysam: \$(pip freeze | grep '^pysam' | sed 's/pysam==//g')
    END_VERSIONS
    """
}
