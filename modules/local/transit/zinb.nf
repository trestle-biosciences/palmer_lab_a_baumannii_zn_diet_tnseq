process TRANSIT_ZINB {
    tag "${contig}_${design_id}"
    label 'process_medium'

    container = 'trestlebiosciences/transit:1.0'

    input:
    path all_design_files
    tuple val(contig), val(design_id), val(design_file), val(ref_condition), path(combined_wig), path(prot), val(covariates)
    path all_wigs    

    output:
    path "*zinb_out.tsv", emit: zinb
    // path "zinb.log", emit: log
    path "versions.yml"           , emit: versions
    path "*.txt", emit: the_wigs

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ? "${contig}_${design_id}-${task.ext.prefix}" : "${contig}_${design_id}"
    def covariates = "${covariates}" ? "--covars ${covariates}" : "" 

    """
    ls *.wig > wigs_in_directory.txt 

    transit zinb \\
        $combined_wig \\
        $design_file \\
        $prot \\
        ${prefix}-zinb_out.tsv \\
        --ref ${ref_condition} \\
        $covariates \\
        $args
        
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        transit: \$(transit --version | sed 's/Version://g')
    END_VERSIONS
    """
}
