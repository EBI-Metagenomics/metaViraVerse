/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { KRONA_KTIMPORTTEXT               } from '../../modules/nf-core/krona/ktimporttext'
include { MULTIQC                          } from '../../modules/nf-core/multiqc/main'

include { EXTRACT_REPS_STATS               } from '../../modules/local/extract_reps_stats'
include { SANKEY_PLOT                      } from '../../modules/local/sankey_plot'

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
        CLUSTERING.out.clusters_tsv,
        gff_and_mapping_all_seqs.map {id, gff, mapfile -> [gff, mapfile]},
    )
    ch_versions = ch_versions.mix(EXTRACT_REPS_STATS.out.versions)

    //
    // Taxonomy visualisation
    //
    KRONA_KTIMPORTTEXT (
        EXTRACT_REPS_STATS.out.reps_krona_tsv
    )
    ch_versions = ch_versions.mix(KRONA_KTIMPORTTEXT.out.versions)

    SANKEY_PLOT (
        EXTRACT_REPS_STATS.out.reps_krona_tsv
    )
    ch_versions = ch_versions.mix(SANKEY_PLOT.out.versions)

    emit:

    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}
