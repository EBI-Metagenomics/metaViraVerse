/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                          } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap                 } from 'plugin/nf-schema'
include { paramsSummaryMultiqc             } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML           } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText           } from '../subworkflows/local/utils_nfcore_metaviraverse_pipeline'

include { PREPROCESSING                    } from '../subworkflows/local/preprocessing'
include { CLUSTERING as CLUSTER_VIRAL_SEQS } from '../subworkflows/local/clustering'
include { CLUSTERING as CLUSTER_PLASMIDS   } from '../subworkflows/local/clustering'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow METAVIRAVERSE {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // Find viral sequences and plasmids
    //
    PREPROCESSING(
       ch_samplesheet
    )

    //
    // Cluster viral sequences
    //
    CLUSTER_VIRAL_SEQS(
       PREPROCESSING.out.viral_seqs,
       95,
       85
    )

    //
    // Cluster plasmids
    //
    CLUSTER_PLASMIDS(
       PREPROCESSING.out.plasmids,
       80,
       85
    )

    //
    // Taxonomy for viral_sequences
    //
    grep from all_gff

    // TODO add VITAP for comparision
    // TODO krona for taxonomy
    // TODO table for reps
    // TODO plots for all reps stats


    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'metaviraverse_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    emit:
    multiqc_report = ''                          //MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
