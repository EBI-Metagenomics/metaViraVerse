/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { KRONA_KTIMPORTTEXT               } from '../../modules/nf-core/krona/ktimporttext'
include { SEQTK_SUBSEQ                     } from '../../modules/nf-core/seqtk/subseq'
include { GUNZIP                           } from '../../modules/nf-core/gunzip'

include { CRISPRCAS_FINDER                 } from '../../modules/local/crispcasfinder'
include { EXTRACT_REPS_STATS               } from '../../modules/local/extract_reps_stats'
include { SANKEY_PLOT                      } from '../../modules/local/sankey_plot'
include { SANKEY_PLOT as SANKEY_VITAP      } from '../../modules/local/sankey_plot'
include { VITAP                            } from '../../modules/local/vitap'

include { CLUSTERING                       } from './clustering'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PROCESS_VIRAL_SEQUENCES {

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
        gff_and_mapping_all_seqs.map {id, gff, mapfile -> [gff, mapfile]}
    )
    ch_versions = ch_versions.mix(EXTRACT_REPS_STATS.out.versions)


    // Extract sequences for cluster reps
    SEQTK_SUBSEQ (
        sequences,
        CLUSTERING.out.clusters_tsv.map{ id, tsv -> tsv }
    )

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

    //
    // Host assignment
    //
    GUNZIP(
        SEQTK_SUBSEQ.out.sequences
    )
    ch_versions = ch_versions.mix(GUNZIP.out.versions)

    CRISPRCAS_FINDER(
        GUNZIP.out.gunzip
    )
    //ch_versions = ch_versions.mix(CRISPRCAS_FINDER.out.versions)

    //
    // Taxonomy VITAP testing...
    //
    VITAP (
        SEQTK_SUBSEQ.out.sequences,
        params.vitap_db
    )
    ch_versions = ch_versions.mix(VITAP.out.versions)

    SANKEY_VITAP (
       VITAP.out.best_lineages
    )
    ch_versions = ch_versions.mix(SANKEY_VITAP.out.versions)

    emit:

    reps_seqs      = SEQTK_SUBSEQ.out.sequences  // compressed
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}
