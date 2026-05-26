/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FIND_CONCATENATE as CONCATENATE_IPHOP_GENOME     } from '../../modules/nf-core/find/concatenate'
include { FIND_CONCATENATE as CONCATENATE_IPHOP_GENUS      } from '../../modules/nf-core/find/concatenate'
include { IPHOP_PREDICT                                    } from '../../modules/nf-core/iphop/predict/main'
include { SEQKIT_SPLIT2 as CHUNK_FNA_IPHOP                 } from '../../modules/nf-core/seqkit/split2'

include { COLLECT_HOST_INFO                                } from '../../modules/local/spacepharer/collect_crispr_host_info'
include { SPACEPHARER_CREATEDB as SPACEPHARER_CREATEDB     } from '../../modules/local/spacepharer/createdb'
include { SPACEPHARER_CREATEDB as SPACEPHARER_CREATEDB_REV } from '../../modules/local/spacepharer/createdb'
include { SPACEPHARER_EASYPREDICT                          } from '../../modules/local/spacepharer/easy_predict'
include { CHANGE_SPACE_TO_UNDERSCORE                       } from '../../modules/local/utils'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow HOST_DETECTION {

    take:
    fna
    metadata

    main:

    ch_versions = channel.empty()
    iphop_host_genome = channel.empty()
    iphop_host_genus = channel.empty()

    if ( !params.skip_iphop) {
        CHUNK_FNA_IPHOP (
            fna,
            [],                                        // length: (disabled) max number of nucleotides per chunk
            params.nucleotide_fasta_chunksize_iphop,   // size: max number of sequences per chunk
        )
        ch_versions = ch_versions.mix(CHUNK_FNA_IPHOP.out.versions)
        def ch_fna_chunks = CHUNK_FNA_IPHOP.out.chunked_output.transpose()

        IPHOP_PREDICT (
            ch_fna_chunks,
            params.iphop_db
        )
        ch_versions = ch_versions.mix(IPHOP_PREDICT.out.versions)

        CONCATENATE_IPHOP_GENOME (
            IPHOP_PREDICT.out.iphop_genome,
            1
        )

        CONCATENATE_IPHOP_GENUS (
            IPHOP_PREDICT.out.iphop_genus,
            1
        )

        iphop_host_genome   = CONCATENATE_IPHOP_GENOME.out.file_out
        iphop_host_genus    = CONCATENATE_IPHOP_GENUS.out.file_out
    }

    if (params.predict_host_from_custom_spacers) {
        CHANGE_SPACE_TO_UNDERSCORE(fna)
        ch_versions = ch_versions.mix(CHANGE_SPACE_TO_UNDERSCORE.out.versions)

        SPACEPHARER_CREATEDB (
            CHANGE_SPACE_TO_UNDERSCORE.out.sanitised_fasta,
            "targetSetDB"
        )
        ch_versions = ch_versions.mix(SPACEPHARER_CREATEDB.out.versions)

        SPACEPHARER_CREATEDB_REV (
            CHANGE_SPACE_TO_UNDERSCORE.out.sanitised_fasta,
            "targetSetDB_rev"
        )
        ch_versions = ch_versions.mix(SPACEPHARER_CREATEDB_REV.out.versions)

        SPACEPHARER_EASYPREDICT(
            channel.of(params.custom_spacers_fasta).map{fasta -> [[id:'spacers'], fasta]},
            "targetSetDB",
            SPACEPHARER_CREATEDB.out.db.map { meta, db -> db },
            SPACEPHARER_CREATEDB_REV.out.db.map { meta, db -> db }
        )
        ch_versions = ch_versions.mix(SPACEPHARER_EASYPREDICT.out.versions)

        COLLECT_HOST_INFO (
           SPACEPHARER_EASYPREDICT.out.predictions,
           channel.of(params.custom_spacers_metadata).map{tsv -> [[id:'spacers'], tsv]},
           metadata
        )
        ch_versions = ch_versions.mix(COLLECT_HOST_INFO.out.versions)
    }

    emit:
    iphop_host_genome   = iphop_host_genome
    iphop_host_genus    = iphop_host_genus
    versions            = ch_versions                 // channel: [ path(versions.yml) ]
}
