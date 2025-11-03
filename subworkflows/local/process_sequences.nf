/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { KRONA_KTIMPORTTEXT               } from '../../modules/nf-core/krona/ktimporttext'
include { MULTIQC                          } from '../../modules/nf-core/multiqc/main'

include { EXTRACT_REPS_STATS               } from '../../modules/local/extract_reps_stats'

include { CLUSTERING                       } from './clustering'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PROCESS_SEQUENCES {

    take:
    sequences
    ani_limit
    coverage_limit
    gff_and_mapping_all_seqs

    main:

    ch_versions = Channel.empty()

    //
    // Cluster sequences
    //
    sequences
        .filter { meta, seqs ->
            seqs && seqs.size() > 0
        }
        .set { samples_seqs }

    CLUSTERING(
       samples_seqs,
       ani_limit,
       coverage_limit
    )
    ch_versions = ch_versions.mix(CLUSTERING.out.versions)

    //
    // Statistics and taxonomy from GFF for viral_sequences reps
    //
    EXTRACT_REPS_STATS (
        gff_and_mapping_all_seqs,
        CLUSTERING.out.clusters_tsv
    )
    ch_versions = ch_versions.mix(EXTRACT_REPS_STATS.out.versions)

    //
    // Taxonomy visualisation
    //
    KRONA_KTIMPORTTEXT (
        EXTRACT_REPS_STATS.out.reps_krona_tsv
    )
    ch_versions = ch_versions.mix(KRONA_KTIMPORTTEXT.out.versions)


    //
    // Taxonomy for viral_sequences
    //
    //grep from all_gff
    // TODO add VITAP for comparision
    // TODO krona for taxonomy
    // TODO table for reps
    // TODO plots for all reps stats

    emit:

    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}
