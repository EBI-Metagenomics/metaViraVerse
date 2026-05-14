/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FIND_CONCATENATE as CONCATENATE_IPHOP_GENOME } from '../../modules/nf-core/find/concatenate'
include { FIND_CONCATENATE as CONCATENATE_IPHOP_GENUS  } from '../../modules/nf-core/find/concatenate'
include { IPHOP_PREDICT                                } from '../../modules/nf-core/iphop/predict/main'
include { SEQKIT_SPLIT2 as CHUNK_FNA_IPHOP             } from '../../modules/nf-core/seqkit/split2'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow HOST_DETECTION {

    take:
    fna

    main:

    ch_versions = channel.empty()

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

    emit:
    iphop_host_genome   = CONCATENATE_IPHOP_GENOME.out.file_out
    iphop_host_genus    = CONCATENATE_IPHOP_GENUS.out.file_out
    versions            = ch_versions                 // channel: [ path(versions.yml) ]
}