process DEEPTOOLS_PLOTPCA {
    tag 'One Process'
    label 'process_low'
    conda (params.enable_conda ? "bioconda::deeptools=3.5.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/deeptools:3.5.1--py_0':
        'quay.io/biocontainers/deeptools:3.5.1--py_0' }"

    input:
    path npz

    output:
    path "*.png", emit: pca_png
    path "*.tsv", emit: pca_tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    
    """
    plotPCA \\
        --corData $npz \\
        --plotFile plotpca.png \
        --outFileNameData plotpca.tsv \\
        $args 
        

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$(plotPCA --version | sed -e "s/plotPCA //g")
    END_VERSIONS
    """
}
