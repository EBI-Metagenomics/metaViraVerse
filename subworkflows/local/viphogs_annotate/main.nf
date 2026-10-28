include { HMMER_HMMSEARCH as HMMER_VIPHOGS   } from '../../../modules/nf-core/hmmer/hmmsearch'
include { FIND_CONCATENATE                   } from '../../../modules/nf-core/find/concatenate'
include { SEQKIT_SPLIT2                      } from '../../../modules/nf-core/seqkit/split2/main'

include { ANNOTATION                         } from '../../../modules/local/annotation'
include { ASSIGN                             } from '../../../modules/local/assign'
include { HMM_POSTPROCESSING                 } from '../../../modules/local/hmm_postprocessing'
include { RATIO_EVALUE                       } from '../../../modules/local/ratio_evalue'

workflow VIPHOGS_ANNOTATION {

    take:
    proteins_faa
    proteins_gff
    viphog_db
    additional_viphog_model_data
    ncbi_db
    viphogs_factor_file

    main:
    // chunk big fasta file
    SEQKIT_SPLIT2(
        proteins_faa,
        [],                                        // length: (disabled) max number of nucleotides per chunk
        params.protein_annotation_fasta_chunksize, // size: max number of sequences per chunk
    )
    def ch_protein_chunks = SEQKIT_SPLIT2.out.chunked_output.transpose()

    hmmer_input = ch_protein_chunks.map{ meta, proteins -> tuple(meta, viphog_db, proteins, false, false, true) }
    HMMER_VIPHOGS(
        hmmer_input
    )

    FIND_CONCATENATE(
        HMMER_VIPHOGS.out.domain_summary.groupTuple(),
        3
    )

    HMM_POSTPROCESSING(
       FIND_CONCATENATE.out.file_out.map { meta, tbl -> record(meta: meta, hmmer_tbl: tbl) }
    )

    // calculate hit qual per protein
    RATIO_EVALUE(
       HMM_POSTPROCESSING.out.map { r -> record(meta: r.meta, modified_table: r.modified_tsv) },
       additional_viphog_model_data,
       params.viphogs_evalue
    )

    // annotate contigs based on ViPhOGs
    ANNOTATION(
       RATIO_EVALUE.out
           .map { r -> [r.meta, r.informative_hits_tsv] }
           .join(proteins_gff)
           .map { meta, tab, gff -> record(meta: meta, tab: tab, gff: gff) }
    )

    // assign lineages
    ASSIGN(
       ANNOTATION.out.map { r -> record(meta: r.meta, tab: r.annotations) },
       ncbi_db,
       viphogs_factor_file
    )

    emit:
    assignment     = ASSIGN.out.map { r -> [r.meta, r.taxonomy_tsv] }
}
