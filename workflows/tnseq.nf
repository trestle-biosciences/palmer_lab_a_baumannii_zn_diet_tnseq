/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters TO-DO: MAKE THIS WORK
// WorkflowTnseq.initialise(params, log)

print params.demultiplex

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
if ( params.demultiplex) {
    def checkPathParamList = [ params.file_to_demultiplex, params.demultiplex_barcodes, 
                               params.multiqc_config, params.fastas_and_gffs ]
    for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }
    } else {
    def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
    for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }
}

// Check mandatory parameters
if ( !params.demultiplex ) {
    if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/* 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
def create_condition_compare_ch(LinkedHashMap row) {
// create meta map
    def comparison = [:]
    comparison.ref         = row.ref
    comparison.condition = row.condition
    return comparison
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { CUTADAPT_DEMULTIPLEX        } from '../modules/local/cutadapt'
include { SEQTK_SAMPLE as SEQTK_RAW   } from '../modules/nf-core/seqtk/sample/main'
include { SEQTK_SAMPLE as SEQTK_TRIM  } from '../modules/nf-core/seqtk/sample/main'
include { FASTQC as FASTQC_RAW        } from '../modules/nf-core/fastqc/main'
include { TRIMGALORE                  } from '../modules/nf-core/trimgalore/main'
include { CUTADAPT                    } from '../modules/nf-core/cutadapt/main'
include { FASTQC as FASTQC_TRIMMED    } from '../modules/nf-core/fastqc/main'
include { CATFILES                    } from '../modules/local/catfiles'
include { CUSTOM_TRANSITPREP          } from '../modules/local/custom/transitprep'
include { GFFREAD                     } from '../modules/nf-core/gffread/main'
include { BOWTIE_ALIGN                } from '../modules/nf-core/bowtie/align/main'
include { BOWTIE_BUILD                } from '../modules/nf-core/bowtie/build/main'
include { SAMTOOLS_SORT               } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX              } from '../modules/nf-core/samtools/index/main'
include { CUSTOM_BAMTOWIG             } from '../modules/local/custom/bamtowig'
include { TRANSIT_CONVERT             } from '../modules/local/transit/convert'
include { TRANSIT_EXPORT              } from '../modules/local/transit/export'
include { TRANSIT_ZINB                } from '../modules/local/transit/zinb'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow TNSEQ {

    ch_versions = Channel.empty()

 
    if (params.demultiplex) {
        
        CUTADAPT_DEMULTIPLEX (
            params.file_to_demultiplex,
            params.demultiplex_barcodes
        )
        ch_versions = ch_versions.mix(CUTADAPT_DEMULTIPLEX.out.versions)

        CUTADAPT_DEMULTIPLEX.out.demultiplexed_fastqs
        .flatMap()
        .map { it -> [['id':"${it}".split("demulti-")[1].replace('.fastq.gz', ''), 
                       'single_end': true] , it] }
        .branch {
        keep: it[0].id !='unknown'
        discard: it[0].id == 'unknown'
        }
        .keep
        .set { ch_reads}       
       
        // CUTADAPT_DEMULTIPLEX.out.demultiplexed_fastqs.map {it -> [meta:it, reads:it] }.view()
    } else {
    //

        // SUBWORKFLOW: Read in samplesheet, validate and stage input files
        //
        INPUT_CHECK (
            ch_input
        )
        ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
        INPUT_CHECK.out.reads.set { ch_reads }
    }

    // CONVERT FASTAS_AND_GFFS PARAMETER TO REQUIRED CHANNELS
        // Read in ids from --input file
    ch_fastas_and_gffs = file(params.fastas_and_gffs)
    ch_designs_file = file(params.designs)
    ch_tests_to_run = file(params.tests_to_run)


    Channel
        .from(ch_fastas_and_gffs)
        .splitCsv(header:true, sep:',', strip:true)
        .map { it -> [it.name, it.fasta, it.gff]}
        .set { ch_references }


    Channel
        .from(ch_fastas_and_gffs)
        .splitCsv(header:true, sep:',', strip:true)
        .map { it -> [it.name, it.fasta]}
        .set { ch_fastas }


    Channel
        .from(ch_fastas_and_gffs)
        .splitCsv(header:true, sep:',', strip:true)
        .map { it -> [it.name, it.gff]}
        .set { ch_gffs }

    ch_fastas
        .map {it -> file(it[1])}
        .collect()
        .set {ch_all_fastas }


    //
    // MODULE: CONCATENATE FASTA FILES 
    //
    CATFILES (
        ch_all_fastas
    )
    //
    // MODULE: Run GFFREAD
    //
    // GFFREAD (
    //     params.gff
    // )
    // ch_versions = ch_versions.mix(GFFREAD.out.versions.first())

    //
    // MODULE: Raw Data Samples
    //
    SEQTK_RAW (
        ch_reads,
        params.sample_size
    )
    ch_versions = ch_versions.mix(SEQTK_RAW.out.versions.first())

    //
    // MODULE: Run FastQC
    //
    FASTQC_RAW (
        ch_reads
    )
    ch_versions = ch_versions.mix(FASTQC_RAW.out.versions.first())


    //
    // MODULE: Run Trimgalore
    //
    // TRIMGALORE (
    //     ch_reads
    // )
    // ch_versions = ch_versions.mix(TRIMGALORE.out.versions)
    //
    // MODULE: Run Cutadapt
    //
    CUTADAPT (
        ch_reads
    )
    ch_versions = ch_versions.mix(CUTADAPT.out.versions.first())

    //
    // MODULE: Run FastQC
    //
    FASTQC_TRIMMED (
        CUTADAPT.out.reads
    )
    ch_versions = ch_versions.mix(FASTQC_TRIMMED.out.versions.first())


    //
    // MODULE: Raw Data Samples
    //
    SEQTK_TRIM (
        CUTADAPT.out.reads,
        params.sample_size
    )
    ch_versions = ch_versions.mix(SEQTK_TRIM.out.versions.first())

    //
    // MODULE: Run Bowtie Build
    //
    BOWTIE_BUILD (
        CATFILES.out.fasta
    )
    ch_versions = ch_versions.mix(BOWTIE_BUILD.out.versions.first())

    //
    // MODULE: Run Bowtie2 Align
    //
    BOWTIE_ALIGN (
        CUTADAPT.out.reads,
        BOWTIE_BUILD.out.index
    )
    ch_versions = ch_versions.mix(BOWTIE_ALIGN.out.versions.first())

    //
    // MODULE: Run samtools sort
    //
    SAMTOOLS_SORT ( 
        BOWTIE_ALIGN.out.bam 
        )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())

    //
    // MODULE: Run samtools index
    //
    SAMTOOLS_INDEX ( 
        SAMTOOLS_SORT.out.bam 
        )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    SAMTOOLS_SORT.out.bam
        .join(SAMTOOLS_INDEX.out.bai, by: [0], remainder: true)
        .join(SAMTOOLS_INDEX.out.csi, by: [0], remainder: true)
        .map {
            meta, bam, bai, csi ->
                if (bai) {
                    [ meta, bam, bai ]
                } else {
                    [ meta, bam, csi ]
                }
        }
        .set { ch_bam_bai }

    //
    // MODULE: Run Deeptools Bamcoverage
    //
    // need to run for each fasta file
    
   
    ch_bam_bai
        .combine(ch_fastas)
        .set{ ch_bam_bai_fasta }

    // ch_bam_bai_fasta.view()


    CUSTOM_BAMTOWIG (
        ch_bam_bai_fasta,
        params.transposon_type,
        params.ta_position
    )
    ch_versions = ch_versions.mix(CUSTOM_BAMTOWIG.out.versions.first())

    CUSTOM_BAMTOWIG.out.transit_wig
        .map {it -> [it[0], it[2], "${file(it[2]).baseName}.wig"]}
        .groupTuple()   
        .map {it -> [it[0], it[1], it[2].join(',')]}
        .set { ch_all_wigs_by_contig }

    CUSTOM_BAMTOWIG.out.transit_wig.collect {it[2] }.set { ch_all_wigs }
    // ch_all_wigs.view()


        // .collect { it[1] }

    // TRANSIT
    // 1. Create Combined Wig File
    // 2. GFF to Prot
    // 3. Run transit zinb 

    TRANSIT_CONVERT (
        ch_gffs
    )
    ch_versions = ch_versions.mix(TRANSIT_CONVERT.out.versions.first())

    TRANSIT_CONVERT
        .out
        .prot
        .map {it -> [it[0], it[1]]}
        .set { ch_prot }

    // Channel.fromPath( params.tests_to_run )
    //     .splitCsv(header: true)
    //     .map { row -> ['design_id': row.design_id, 
    //                    'ref_condition': "${row.ref_condition}",
    //                    'samples': "${row.samples}".replaceAll('-', ',')
    //     ]}
    //     .set {ch_designs}

    //                     // 'wigs': "${row.wigs}".replaceAll('-', ','),
    //                     // 'covariates': row.covariates ?: ''

    // ch_designs
    //     .combine(ch_prot)
    //     .map { it ->   [it[1], //contig name
    //                     it[0].design_id,
    //                     it[0].ref_condition,
    //                     it[0].samples, 
    //                     file(it[2]) //prot file
    //                     ]} 
    //     .set { ch_designs }
  

    // contig name, design_id, ref_condition, prot file, all wig files, wig files as a list
    // ch_designs.join(ch_all_wigs)
    //     .set {ch_transit_prep }

    CUSTOM_TRANSITPREP (
        ch_fastas_and_gffs,
        ch_designs_file,
        ch_tests_to_run
    )
    ch_versions = ch_versions.mix(CUSTOM_TRANSITPREP.out.versions.first())

    CUSTOM_TRANSITPREP.out
        .all_designs
        .splitCsv(header: true)
        .map { row -> [row.contig,
                       row.design_id,
                       row.design_file,
                       row.ref_condition,
                       "${row.wigs}".replaceAll('-', ','),
                       row.covariates
        ]}
        .set { ch_designs }   

    ch_prot
        .cross( ch_designs )
        .map { it -> [it[1][0], // CONTIG NAME
                      it[1][1], // DESIGN ID
                      it[1][2], // DESIGN FILE
                      it[1][3], // REFERENCE CONDITION
                      it[1][4], // WIG FILES
                      it[0][1], // PROT FILE
                      it[1][5]  // COVARIATES
        ]}
        .set { ch_designs }
        
    TRANSIT_EXPORT (
        ch_all_wigs,
        ch_designs
    )
    ch_versions = ch_versions.mix(TRANSIT_EXPORT.out.versions.first())

    // TRANSIT_EXPORT.out.combined_wigs.view()

    // CREATE CHANNEL OF COMPARISONS TO RUN
    // def conditions_to_compare_file = file(params.conditions_to_compare)
    // ch_cc_file = Channel.fromPath(conditions_to_compare_file)


    
    // print(conditions_to_compare_file.splitCsv ( header:true, sep:'\t' ))

    // ch_cc_file
    //     .splitCsv ( header:true, sep:'\t' )
    //     .map { it -> [it.ref,it.include_conditions, it.exclude_conditions] }
    //     .set { ch_conditions_to_compare }

    // ch_conditions_to_compare.view()

 
    TRANSIT_ZINB (
        CUSTOM_TRANSITPREP.out.all_design_files,
        TRANSIT_EXPORT.out.combined_wigs,
        ch_all_wigs
    )
    ch_versions = ch_versions.mix(TRANSIT_ZINB.out.versions.first())

    //
    // MODULE: CUSTOM_DUMPSOFTWAREVERSIONS
    //
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowTnseq.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    // ch_multiqc_files = ch_multiqc_files.mix(FASTQC_RAW.out.zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_TRIMMED.out.zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(CUTADAPT.out.log.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(BOWTIE_ALIGN.out.log.collect{it[1]}.ifEmpty([]))
    // ch_multiqc_files = ch_multiqc_files.mix(TRIMGALORE.out.log.collect{it[1]}.ifEmpty([]))
    // ch_multiqc_files = ch_multiqc_files.mix(TRIMGALORE.out.zip.collect{it[1]}.ifEmpty([]))
    // ch_multiqc_files = ch_multiqc_files.mix(TRANSIT_ZINB.out.log.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect()
    )
    multiqc_report = MULTIQC.out.report.toList()
    ch_versions    = ch_versions.mix(MULTIQC.out.versions)


}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
