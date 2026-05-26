include { HMMER_HMMSEARCH as HMMER_VIPHOGS   } from '../../modules/nf-core/hmmer/hmmsearch'
include { FIND_CONCATENATE                   } from '../../modules/nf-core/find/concatenate'
include { SEQKIT_SPLIT2                      } from '../../modules/nf-core/seqkit/split2/main'

include { ANNOTATION                         } from '../../modules/local/annotation'
include { ASSIGN                             } from '../../modules/local/assign'
include { HMM_POSTPROCESSING                 } from '../../modules/local/hmm_postprocessing'
include { RATIO_EVALUE                       } from '../../modules/local/ratio_evalue'

workflow VIPHOGS_ANNOTATION {

    take:
    proteins_faa
    proteins_gff
    viphog_db
    additional_model_data
    ncbi_db
    factor_file

    main:
    // chunk big fasta file
    SEQKIT_SPLIT2(
        proteins_faa,
        [],                                        // length: (disabled) max number of nucleotides per chunk
        params.protein_annotation_fasta_chunksize, // size: max number of sequences per chunk
    )
    def ch_protein_chunks = SEQKIT_SPLIT2.out.chunked_output.transpose()

    hmmer_input = ch_protein_chunks.map{ meta, proteins -> tuple(meta, viphog_db, proteins, false, true, false) }
    HMMER_VIPHOGS(
        hmmer_input
    )

    FIND_CONCATENATE(
        HMMER_VIPHOGS.out.target_summary.groupTuple(by: [0,1]),
        3
    )

    HMM_POSTPROCESSING(
       FIND_CONCATENATE.out.file_out
    )

    // calculate hit qual per protein
    RATIO_EVALUE(
       HMM_POSTPROCESSING.out,
       additional_model_data
    )

    // annotate contigs based on ViPhOGs
    ANNOTATION(
       RATIO_EVALUE.out.join(proteins_gff, by:[0,1])
    )

    // assign lineages
    factor_file = Channel.empty()
    if (params.factor_file) {
        factor_file = file(params.factor_file, checkIfExists: true)
    }
    ASSIGN(
       ANNOTATION.out.annotations,
       ncbi_db,
       factor_file
    )
}