/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { SEQTK_SUBSEQ                     } from '../../modules/nf-core/seqtk/subseq'

include { CLUSTERING                       } from './clustering'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PROCESS_PLASMIDS {

    take:
    sequences

    main:

    ch_versions = channel.empty()

    //
    // Cluster sequences
    //
    sequences
        .filter { _meta, seqs ->
            seqs && seqs.size() > 0
        }
        .set { samples_seqs }

    CLUSTERING(
       samples_seqs,
       'gani',
       false,
       0.35,
       false
    )
    ch_versions = ch_versions.mix(CLUSTERING.out.versions)

    if ( params.cluster_vclust ) {
        // Take second column unique values
        CLUSTERING.out.clusters_tsv
            .splitCsv(sep: '\t', skip: 0)  // Adjust separator and skip header
            .map { row -> row[1] }          // Extract second column (0-indexed)
            .unique()                        // Get unique values
            .collectFile(name: 'plasmid_reps_vclust.txt', newLine: true)  // Write to file
            .set { reps_ch }
    } else {
        reps_ch = CLUSTERING.out.clusters_tsv.map{ _id, tsv -> tsv }  // for blastn all reps are in first column
    }
    // Extract sequences for cluster reps
    SEQTK_SUBSEQ (
        sequences,
        reps_ch
    )
    ch_versions = ch_versions.mix(SEQTK_SUBSEQ.out.versions)

    emit:

    reps_tsv       = reps_ch
    reps_seqs      = SEQTK_SUBSEQ.out.sequences  // compressed
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}
