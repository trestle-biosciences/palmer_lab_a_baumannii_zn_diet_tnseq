process TRANSIT_EXPORT {
    tag "${contig}_${design_id}"
    label 'process_low'

    // container = 'bradfordtw/transit:1.0'

    conda (params.enable_conda ? "bioconda::transit=3.2.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/transit:3.2.3--pyhdfd78af_0':
        'quay.io/biocontainers/transit:3.2.3--pyhdfd78af_0' }"

    input:
    path all_wigs    
    tuple val(contig), val(design_id), val(design_file), val(ref_condition), val(wigs), path(prot), val(covariates)

    output:
    tuple val(contig), val(design_id), val(design_file), val(ref_condition), path("*.combined_wig.wig"), path(prot), val(covariates), emit: combined_wigs
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    
    """
    transit \\
        export \\
        combined_wig \\
        $wigs \\
        $prot \\
        ${contig}_${design_id}.combined_wig.wig \\
        $args 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        transit: \$(transit --version | sed 's/Version://g')
    END_VERSIONS
    """
}
