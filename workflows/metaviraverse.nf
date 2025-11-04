/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap                      } from 'plugin/nf-schema'
include { paramsSummaryMultiqc                  } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                } from '../subworkflows/local/utils_nfcore_metaviraverse_pipeline'

include { PREPROCESSING                         } from '../subworkflows/local/preprocessing'
include { PROCESS_SEQUENCES as PROCESS_VIRUSES  } from '../subworkflows/local/process_sequences'
include { PROCESS_SEQUENCES as PROCESS_PLASMIDS } from '../subworkflows/local/process_sequences'

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
    ch_versions = ch_versions.mix(PREPROCESSING.out.versions)

    //
    // Process viral sequences
    //
    PROCESS_VIRUSES(
       PREPROCESSING.out.viral_seqs,
       95,
       85,
       PREPROCESSING.out.all_gff.join( PREPROCESSING.out.all_mapping )
    )
    ch_versions = ch_versions.mix(PROCESS_VIRUSES.out.versions)

    //
    // Process plasmids
    //
    PROCESS_PLASMIDS(
       PREPROCESSING.out.plasmids,
       80,
       85,
       PREPROCESSING.out.all_gff.join( PREPROCESSING.out.all_mapping )
    )
    ch_versions = ch_versions.mix(PROCESS_PLASMIDS.out.versions)

    //
    // Taxonomy for viral_sequences
    //
    //grep from all_gff
    // TODO at all gffs
    // TODO cat all mapfiles
    // TODO add VITAP for comparision
    // TODO sankey for taxonomy



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
