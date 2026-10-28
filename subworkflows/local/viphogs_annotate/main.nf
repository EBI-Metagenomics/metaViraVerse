nextflow.enable.types = true

include { HMMER_HMMSEARCH as HMMER_VIPHOGS   } from '../../../modules/nf-core/hmmer/hmmsearch'
include { FIND_CONCATENATE                   } from '../../../modules/nf-core/find/concatenate'
include { SEQKIT_SPLIT2                      } from '../../../modules/nf-core/seqkit/split2/main'

include { ANNOTATION                         } from '../../../modules/local/annotation'
include { ASSIGN                             } from '../../../modules/local/assign'
include { HMM_POSTPROCESSING                 } from '../../../modules/local/hmm_postprocessing'
include { RATIO_EVALUE                       } from '../../../modules/local/ratio_evalue'

record ViphogsParams {
    viphog_db: Path
    additional_viphog_model_data: Path
    ncbi_db: Path
    viphogs_factor_file: Path?
    viphogs_evalue: Float
    protein_annotation_fasta_chunksize: Integer
}

record ViphogsInput {
    meta: Map
    faa: Path   // representative proteins, uncompressed
    gff: Path   // representative GFF
}

record ViphogsAssignment {
    meta: Map
    taxonomy_tsv: Path
}

workflow VIPHOGS_ANNOTATION {

    take:
    proteins: Channel<ViphogsInput>
    opts: ViphogsParams

    main:
    // chunk big fasta file
    chunks = SEQKIT_SPLIT2(
        proteins.map { r -> tuple(r.meta, r.faa) },
        [],                                        // length: (disabled) max number of nucleotides per chunk
        opts.protein_annotation_fasta_chunksize,   // size: max number of sequences per chunk
    )
    // a single chunk is emitted as a path, several as a list: transpose handles both
    def ch_protein_chunks = chunks.chunked_output.transpose()

    hmmer = HMMER_VIPHOGS(
        ch_protein_chunks.map { meta, faa -> tuple(meta, opts.viphog_db, faa, false, false, true) }
    )

    hmmer_domtbl = FIND_CONCATENATE(
        hmmer.domain_summary.groupTuple(),
        3
    )

    postprocessed = HMM_POSTPROCESSING(
       hmmer_domtbl.file_out.map { meta, tbl -> record(meta: meta, hmmer_tbl: tbl) }
    )

    // calculate hit qual per protein
    ratio = RATIO_EVALUE(
       postprocessed.map { r -> record(meta: r.meta, modified_table: r.modified_tsv) },
       opts.additional_viphog_model_data,
       opts.viphogs_evalue
    )

    // annotate contigs based on ViPhOGs
    annotation = ANNOTATION(
       ratio.join(proteins, by: 'meta')
           .map { r -> record(meta: r.meta, tab: r.informative_hits_tsv, gff: r.gff) }
    )

    // assign lineages
    assignment = ASSIGN(
       annotation.map { r -> record(meta: r.meta, tab: r.annotations) },
       opts.ncbi_db,
       opts.viphogs_factor_file
    )

    emit:
    assignment: Channel<ViphogsAssignment> = assignment
}
