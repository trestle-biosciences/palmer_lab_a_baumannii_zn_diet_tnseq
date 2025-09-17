process TRANSIT_CONVERT {
    tag "$gff"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::transit=3.2.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/transit:3.2.3--pyhdfd78af_0':
        'quay.io/biocontainers/transit:3.2.3--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(gff)

    output:
    tuple val(meta), path('*.prot'), emit: prot
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    
    
    """
    transit \\
        convert  \\
        gff_to_prot_table \\
        $gff \\
        ${meta}_prot_for_transit.prot \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        transit: \$(transit --version | sed 's/Version://g')
    END_VERSIONS
    """
}
