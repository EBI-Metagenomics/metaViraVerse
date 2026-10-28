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


record ProteinsParams {
    amrfinderplus_db: String?
    deeparg_db: String?
    deeparg_db_version: String?
    deeparg_model: String?
    deeparg_tool_version: String?
    rgi_db: String?
    mmseqs_id_threshold: Float
    mmseqs_cov_threshold: Float
    protein_annotation_fasta_chunksize: Integer
}

record ProteinsInput {
    meta: Map
    faa: Path   // representative proteins, uncompressed
    gff: Path   // representative GFF
}

record AmrGff {
    meta: Map
    amr_gff: Path
}

record AmrReport {
    meta: Map
    report: Path
}

record HmmerTables {
    meta: Map
    hmmer_tables: List<Path>
    hmmer_summaries: List<Path>
}

record ClusterPlots {
    meta: Map
    size_distribution_tsv: Path
    size_distribution_png: Path
    rank_abundance: Path
    top_functions_tsv: Path
    top_functions_png: Path?
    purity_tsv: Path
    purity_png: Path?
    stats: Path
}

record Phams {
    meta: Map
    phams_tsv: Path
    phams_plots: ClusterPlots
}

record MmseqsClusters {
    meta: Map
    mmseqs_tsv: Path
    mmseqs_plots: ClusterPlots
}

// Optional steps leave their fields null
record ProteinsResult {
    meta: Map
    hmmer_tables: List<Path>?      // null when --annotation_db is not set
    hmmer_summaries: List<Path>?   // per-database function summaries
    amr_gff: Path?                 // null when all AMR tools are skipped
    amrfinderplus_tsv: Path?
    deeparg_tsv: Path?             // hAMRonization report
    rgi_tsv: Path?                 // hAMRonization report
    phams_tsv: Path?
    phams_plots: ClusterPlots?
    mmseqs_tsv: Path?
    mmseqs_plots: ClusterPlots?
}

workflow PROTEINS_PROCESSING {

    take:
    input: Channel<ProteinsInput>
    skip_amrfinderplus: Boolean
    skip_deeparg: Boolean
    skip_rgi: Boolean
    skip_phammseqs: Boolean
    skip_mmseqs: Boolean
    hmm_db: Path?     // directory with *.hmm.gz and *.tsv; null to skip HMMER annotation
    opts: ProteinsParams

    main:

    def amr_gff: Channel<AmrGff> = channel.empty()
    def amrfinderplus: Channel<AmrReport> = channel.empty()
    def deeparg: Channel<AmrReport> = channel.empty()
    def rgi: Channel<AmrReport> = channel.empty()
    def hmmer_tables: Channel<HmmerTables> = channel.empty()
    def phams: Channel<Phams> = channel.empty()
    def mmseqs_clusters: Channel<MmseqsClusters> = channel.empty()

    //
    // -------- Antimicrobial resistence detection
    //
    if (!(skip_amrfinderplus && skip_deeparg && skip_rgi)) {
        amr = AMR_ANNOTATION (
            input.map { r -> tuple(r.meta, r.faa, r.gff) },
            opts.amrfinderplus_db,
            opts.deeparg_db,
            opts.deeparg_db_version,
            opts.deeparg_model,
            opts.deeparg_tool_version,
            opts.rgi_db,
            skip_amrfinderplus,
            skip_deeparg,
            skip_rgi
        )
        amr_gff = amr.gff.map { meta, gff -> record(meta: meta, amr_gff: gff) }
        // tool reports carry `is_proteins` in their meta; read them only when the tool ran
        // (an empty channel from the untyped subworkflow cannot be read here)
        if ( !skip_amrfinderplus ) {
            amrfinderplus = amr.amrfinderplus.map { meta, f -> record(meta: meta.subMap(meta.keySet().findAll { k -> k != 'is_proteins' }), report: f) }
        }
        if ( !skip_deeparg ) {
            deeparg = amr.deeparg.map { meta, f -> record(meta: meta.subMap(meta.keySet().findAll { k -> k != 'is_proteins' }), report: f) }
        }
        if ( !skip_rgi ) {
            rgi = amr.rgi.map { meta, f -> record(meta: meta.subMap(meta.keySet().findAll { k -> k != 'is_proteins' }), report: f) }
        }
    }

    if ( !skip_phammseqs ) {
        //
        // -------- Assort phage protein sequences into phamilies using MMseqs2
        //
        phammseqs = PHAMMSEQS(
           input.map { r -> record(meta: r.meta, faa: r.faa) }
        )

        phams_plots = PLOT_PHAMS_STATS (
            // functions are already in the phams TSV; no reference labelling: compare phams with protein functions
            phammseqs.map { r -> record(meta: r.meta, clusters: r.phams_tsv, fasta: null, reference: null) },
            'phammseqs'
        )
        phams = phammseqs.join(phams_plots.map { r -> record(meta: r.meta, phams_plots: r) }, by: 'meta')
            .map { r -> record(meta: r.meta, phams_tsv: r.phams_tsv, phams_plots: r.phams_plots) }
    }

    if ( !skip_mmseqs ) {
        mmseqs = MMSEQS (
            input.map { r -> record(meta: r.meta, faa_file: r.faa) },
            opts.mmseqs_id_threshold,
            opts.mmseqs_cov_threshold
        )

        // take functions from the protein FASTA headers
        mmseqs_plots = PLOT_MMSEQS_STATS (
            mmseqs.join(input, by: 'meta')
                .map { r -> record(meta: r.meta, clusters: r.mmseq_cluster_tsv, fasta: r.faa, reference: null) },
            'mmseqs'
        )
        mmseqs_clusters = mmseqs.join(mmseqs_plots.map { r -> record(meta: r.meta, mmseqs_plots: r) }, by: 'meta')
            .map { r -> record(meta: r.meta, mmseqs_tsv: r.mmseq_cluster_tsv, mmseqs_plots: r.mmseqs_plots) }
    }

    if (hmm_db) {
        //
        // -------- Protein chunking for parallel annotation
        // Default: 50,000 sequences per chunk (override with --protein_annotation_fasta_chunksize)
        //

        chunks = SEQKIT_SPLIT2(
            input.map { r -> tuple(r.meta, r.faa) },
            [],                                        // length: (disabled) max number of nucleotides per chunk
            opts.protein_annotation_fasta_chunksize,   // size: max number of sequences per chunk
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

        summaries = SUMMARISE_ANNOTATIONS(
            hmmer_tblout.file_out
                .join(hmm_metadata, by: 0)
                .map { meta, tbl, metadata -> record(meta: meta, tbl: tbl, metadata: metadata) }
        )

        // all tables, under the input's meta (viruses/plasmids)
        def tables = hmmer_tblout.file_out
            .map { _meta, table -> table }
            .collect()
            .map { ts -> ts.toList() }
        def summary_tables = summaries
            .map { r -> r.summary_tsv }
            .collect()
            .map { ts -> ts.toList() }
        hmmer_tables = input
            .map { r -> record(meta: r.meta) }
            .combine(hmmer_tables: tables)
            .combine(hmmer_summaries: summary_tables)
    }

    // Steps above are optional: their fields are null when skipped
    result = input
        .map { r -> record(meta: r.meta) }
        .join(hmmer_tables, by: 'meta', remainder: true)
        .join(amr_gff, by: 'meta', remainder: true)
        .join(amrfinderplus.map { r -> record(meta: r.meta, amrfinderplus_tsv: r.report) }, by: 'meta', remainder: true)
        .join(deeparg.map { r -> record(meta: r.meta, deeparg_tsv: r.report) }, by: 'meta', remainder: true)
        .join(rgi.map { r -> record(meta: r.meta, rgi_tsv: r.report) }, by: 'meta', remainder: true)
        .join(phams, by: 'meta', remainder: true)
        .join(mmseqs_clusters, by: 'meta', remainder: true)
        .map { r ->
            record(
                meta: r.meta,
                hmmer_tables: r.hmmer_tables, hmmer_summaries: r.hmmer_summaries,
                amr_gff: r.amr_gff, amrfinderplus_tsv: r.amrfinderplus_tsv, deeparg_tsv: r.deeparg_tsv, rgi_tsv: r.rgi_tsv,
                phams_tsv: r.phams_tsv, phams_plots: r.phams_plots,
                mmseqs_tsv: r.mmseqs_tsv, mmseqs_plots: r.mmseqs_plots
            )
        }

    emit:
    result: Channel<ProteinsResult> = result
}
