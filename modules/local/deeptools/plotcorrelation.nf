process DEEPTOOLS_PLOTCORRELATION {
    tag 'One Process'
    label 'process_low'
    conda (params.enable_conda ? "bioconda::deeptools=3.5.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/deeptools:3.5.1--py_0':
        'quay.io/biocontainers/deeptools:3.5.1--py_0' }"

    input:
    path npz
    val cor_method

    output:
    path "*.png", emit: cor_png
    path "*.tsv", emit: cor_tsv

    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    plotCorrelation \\
        --corData  $npz \\
        --corMethod $cor_method \\
        --whatToPlot 'heatmap' \\
        --plotFile correlation_heatmap.png \\
        --outFileCorMatrix correlation_heatmap.tsv \\
        $args 


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
