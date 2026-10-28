nextflow.enable.types = true

include { MMSEQS                                                  } from '../../../modules/local/mmseqs'
include { PHAMMSEQS                                               } from '../../../modules/local/phammseqs'
include { SUMMARISE_ANNOTATIONS                                   } from '../../../modules/local/summarise_annotations'
include { PLOT_PROTEIN_CLUSTERS_STATS as PLOT_PHAMS_STATS         } from '../../../modules/local/plot_protein_clusters_stats'
include { PLOT_PROTEIN_CLUSTERS_STATS as PLOT_MMSEQS_STATS        } from '../../../modules/local/plot_protein_clusters_stats'

include { HMMER_HMMSEARCH                                         } from '../../../modules/nf-core/hmmer/hmmsearch'
include { SEQKIT_SPLIT2                                           } from '../../../modules/nf-core/seqkit/split2'
include { FIND_CONCATENATE as CONCATENATE_HMMER_TBLOUT            } from '../../../modules/nf-core/find/concatenate'

include { AMR_ANNOTATION                                          } from '../../ebi-metagenomics/amr_annotation'


record ProteinsInput {
    meta: Map
    faa: Path   // representative proteins, uncompressed
    gff: Path   // representative GFF
}

record AmrGff {
    meta: Map
    amr_gff: Path
}

record HmmerTables {
    meta: Map
    hmmer_tables: List<Path>
}

record ProteinsResult {
    meta: Map
    hmmer_tables: List<Path>?   // null when --annotation_db is not set
    amr_gff: Path?              // null when all AMR tools are skipped
}

workflow PROTEINS_PROCESSING {

    take:
    input: Channel<ProteinsInput>
    skip_amrfinderplus: Boolean
    skip_deeparg: Boolean
    skip_rgi: Boolean
    skip_phammseqs: Boolean
    skip_mmseqs: Boolean
    hmm_db: String?   // directory with *.hmm.gz and *.tsv; null to skip HMMER annotation

    main:

    def amr_gff: Channel<AmrGff> = channel.empty()
    def hmmer_tables: Channel<HmmerTables> = channel.empty()

    //
    // -------- Antimicrobial resistence detection
    //
    if (!(skip_amrfinderplus && skip_deeparg && skip_rgi)) {
        amr = AMR_ANNOTATION (
            input.map { r -> tuple(r.meta, r.faa, r.gff) },
            params.amrfinderplus_db,
            params.deeparg_db,
            params.deeparg_db_version,
            params.deeparg_model,
            params.deeparg_tool_version,
            params.rgi_db,
            skip_amrfinderplus,
            skip_deeparg,
            skip_rgi
        )
        // a single-emit (untyped) workflow returns that channel itself
        amr_gff = amr.map { meta, gff -> record(meta: meta, amr_gff: gff) }
    }

    if ( !skip_phammseqs ) {
        //
        // -------- Assort phage protein sequences into phamilies using MMseqs2
        //
        phams = PHAMMSEQS(
           input.map { r -> record(meta: r.meta, faa: r.faa) }
        )

        PLOT_PHAMS_STATS (
            // functions are already in the phams TSV; no reference labelling: compare phams with protein functions
            phams.map { r -> record(meta: r.meta, clusters: r.phams_tsv, fasta: null, reference: null) },
            'phammseqs'
        )
    }

    if ( !skip_mmseqs ) {
        mmseqs = MMSEQS (
            input.map { r -> record(meta: r.meta, faa_file: r.faa) },
            params.mmseqs_id_threshold,
            params.mmseqs_cov_threshold
        )

        // take functions from the protein FASTA headers
        PLOT_MMSEQS_STATS (
            mmseqs.join(input, by: 'meta')
                .map { r -> record(meta: r.meta, clusters: r.mmseq_cluster_tsv, fasta: r.faa, reference: null) },
            'mmseqs'
        )
    }

    if (hmm_db) {
        //
        // -------- Protein chunking for parallel annotation
        // Default: 50,000 sequences per chunk (override with --protein_annotation_fasta_chunksize)
        //

        chunks = SEQKIT_SPLIT2(
            input.map { r -> tuple(r.meta, r.faa) },
            [],                                        // length: (disabled) max number of nucleotides per chunk
            params.protein_annotation_fasta_chunksize, // size: max number of sequences per chunk
        )

        // a single chunk is emitted as a path, several as a list: transpose handles both
        def ch_protein_chunks = chunks.chunked_output.transpose()

        //
        // -------- Functional annotation with HMMER
        //
        def hmm_files = channel.fromPath("${hmm_db}/*.hmm.gz").collect()

        def hmmsearch_input = ch_protein_chunks
            .map { _meta, faa_file -> record(faa_file: faa_file) }
            .combine(hmm_files: hmm_files)
            .flatMap { r ->
                r.hmm_files.collect { hmm_file ->
                    def meta = [
                        id: hmm_file.baseName.replace('.hmm', '')
                    ]
                    tuple(meta, hmm_file, r.faa_file, false, true, false)
                }
            }

        hmmsearch = HMMER_HMMSEARCH (
            hmmsearch_input
        )

        hmmer_tblout = CONCATENATE_HMMER_TBLOUT (
            hmmsearch.target_summary.groupTuple(),
            3
        )

        //
        // ----- Add metadata and generate a full summary
        //
        def hmm_metadata = channel
            .fromPath("${hmm_db}/*.tsv")
            .map { hmm_file ->
                 def meta = [
                    id: hmm_file.baseName
                 ]
                 tuple(meta, hmm_file)
            }

        SUMMARISE_ANNOTATIONS(
            hmmer_tblout.file_out
                .join(hmm_metadata, by: 0)
                .map { meta, tbl, metadata -> record(meta: meta, tbl: tbl, metadata: metadata) }
        )

        // all tables, under the input's meta (viruses/plasmids)
        def tables = hmmer_tblout.file_out
            .map { _meta, table -> table }
            .collect()
        hmmer_tables = input
            .map { r -> record(meta: r.meta) }
            .combine(hmmer_tables: tables)
    }

    // Steps above are optional: their fields are null when skipped
    result = input
        .map { r -> record(meta: r.meta) }
        .join(hmmer_tables, by: 'meta', remainder: true)
        .join(amr_gff, by: 'meta', remainder: true)
        .map { r -> record(meta: r.meta, hmmer_tables: r.hmmer_tables, amr_gff: r.amr_gff) }

    emit:
    result: Channel<ProteinsResult> = result
}
