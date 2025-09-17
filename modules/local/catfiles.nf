process CATFILES {
    tag 'All Files'
    label 'process_single'

    conda (params.enable_conda ? "bioconda::samtools=1.15.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.15.1--h1170115_0' :
        'quay.io/biocontainers/samtools:1.15.1--h1170115_0' }"

    input:
    path list_of_files

    output:
    path "combined_contigs.fasta", emit: fasta

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    
    """
    cat *.fasta > combined_contigs.fasta

    """
}
   