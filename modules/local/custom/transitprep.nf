process CUSTOM_TRANSITPREP {
    tag "all_designs"
    label 'process_single'
   
    conda (params.enable_conda ? "bioconda pandas==1.4.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.4.3':
        'quay.io/biocontainers/pandas:1.4.3' }"

    input:
    path fastas_and_gffs
    path design_file
    path tests_to_run

    output:
    path ("all_designs.txt"), emit: all_designs
    path ("*design.txt"), emit: all_design_files
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    transit_prep.py \\
        --design-file $design_file \\
        --tests-to-run $tests_to_run \\
        --fastas-and-gffs $fastas_and_gffs \\
        --outfile all_designs.txt \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(pip freeze | grep '^pandas' | sed 's/pandas==//g')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${contig}_${design_id}"
    """
    touch ${prefix}_design.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(python --version))
    END_VERSIONS
    """
}
