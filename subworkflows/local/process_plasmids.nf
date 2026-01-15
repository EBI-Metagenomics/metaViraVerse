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
        .filter { meta, seqs ->
            seqs && seqs.size() > 0
        }
        .set { samples_seqs }

    CLUSTERING(
       samples_seqs,
       params.blastn_ani_threshold_plasmid,
       params.blastn_cov_threshold_plasmid
    )
    ch_versions = ch_versions.mix(CLUSTERING.out.versions)


    // Extract sequences for cluster reps
    SEQTK_SUBSEQ (
        sequences,
        CLUSTERING.out.clusters_tsv.map{ id, tsv -> tsv }
    )

    emit:

    reps_seqs      = SEQTK_SUBSEQ.out.sequences  // compressed
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}
