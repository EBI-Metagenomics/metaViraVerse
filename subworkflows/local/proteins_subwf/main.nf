include { MMSEQS                                                  } from '../../../modules/local/mmseqs'
include { PHAMMSEQS                                               } from '../../../modules/local/phammseqs'
include { SUMMARISE_ANNOTATIONS                                   } from '../../../modules/local/summarise_annotations'
include { PLOT_PROTEIN_CLUSTERS_STATS as PLOT_PHAMS_STATS         } from '../../../modules/local/plot_protein_clusters_stats'
include { PLOT_PROTEIN_CLUSTERS_STATS as PLOT_MMSEQS_STATS        } from '../../../modules/local/plot_protein_clusters_stats'

include { HMMER_HMMSEARCH                                         } from '../../../modules/nf-core/hmmer/hmmsearch'
include { SEQKIT_SPLIT2                                           } from '../../../modules/nf-core/seqkit/split2'
include { FIND_CONCATENATE as CONCATENATE_HMMER_TBLOUT            } from '../../../modules/nf-core/find/concatenate'

include { AMR_ANNOTATION                                          } from '../../ebi-metagenomics/amr_annotation'


workflow PROTEINS_PROCESSING {

    take:
    input  // (meta, reps_faa.uncompressed, reps_gff)
    skip_amrfinderplus
    skip_deeparg
    skip_rgi
    skip_phammseqs
    skip_mmseqs
    hmm_db

    main:

    ch_versions  = channel.empty()
    hmmer_tables = channel.empty()
    amr_gff      = channel.empty()

    ch_proteins  = input.map{ meta, faa, _gff -> tuple(meta, faa) }

    //
    // -------- Antimicrobial resistence detection
    //
    if (!(skip_amrfinderplus && skip_deeparg && skip_rgi)) {
        AMR_ANNOTATION (
            input,
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
        amr_gff = AMR_ANNOTATION.out.gff
    }

    if ( !skip_phammseqs ) {
        //
        // -------- Assort phage protein sequences into phamilies using MMseqs2
        //
        PHAMMSEQS(
           ch_proteins
        )
        ch_versions = ch_versions.mix(PHAMMSEQS.out.versions)

        PLOT_PHAMS_STATS (
            PHAMMSEQS.out.phams_tsv,
            'phammseqs',
            [],   // functions are already in the phams TSV
            []    // no reference labelling: compare phams with protein functions
        )
        ch_versions = ch_versions.mix(PLOT_PHAMS_STATS.out.versions)
    }

    if ( !skip_mmseqs ) {
        MMSEQS (
            ch_proteins,
            channel.value(params.mmseqs_id_threshold),
            channel.value(params.mmseqs_cov_threshold)
        )

        // MMSEQS emits the cluster TSV without meta: attach the proteins' meta,
        // and take functions from the protein FASTA headers
        PLOT_MMSEQS_STATS (
            ch_proteins.map { meta, _faa -> meta }.combine(MMSEQS.out.mmseq_cluster_tsv),
            'mmseqs',
            ch_proteins.map { _meta, faa -> faa },
            []
        )
        ch_versions = ch_versions.mix(PLOT_MMSEQS_STATS.out.versions)
    }

    if (hmm_db) {
        //
        // -------- Protein chunking for parallel annotation
        // Default: 50,000 sequences per chunk (override with --protein_annotation_fasta_chunksize)
        //

        SEQKIT_SPLIT2(
            ch_proteins,
            [],                                        // length: (disabled) max number of nucleotides per chunk
            params.protein_annotation_fasta_chunksize, // size: max number of sequences per chunk
        )
        ch_versions = ch_versions.mix(SEQKIT_SPLIT2.out.versions)

        def ch_protein_chunks = SEQKIT_SPLIT2.out.chunked_output.transpose()

        //
        // -------- Functional annotation with HMMER
        //
        def hmm_ch = channel
            .fromPath("${hmm_db}/*.hmm.gz")
            .map { hmm_file -> tuple(hmm_file.baseName, hmm_file) }

        def hmmsearch_input = hmm_ch
            .combine(ch_protein_chunks)
            .map { hmm_id, hmm_file, faa_id, faa_file ->
                def meta = [
                    id: hmm_id.replace('.hmm', '')
                ]
                tuple(meta, hmm_file, faa_file, false, true, false)
            }

        HMMER_HMMSEARCH (
            hmmsearch_input
        )
        ch_versions = ch_versions.mix(HMMER_HMMSEARCH.out.versions)

        CONCATENATE_HMMER_TBLOUT (
            HMMER_HMMSEARCH.out.target_summary.groupTuple(),
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
            CONCATENATE_HMMER_TBLOUT.out.file_out.join(hmm_metadata)
        )
        ch_versions = ch_versions.mix(SUMMARISE_ANNOTATIONS.out.versions)

        hmmer_tables = CONCATENATE_HMMER_TBLOUT.out.file_out
            .map { _meta, table -> table }                                      // Extract just the table files
            .collect()                                                          // Gather all tables into a list
            .combine(ch_proteins.map { meta, _faa -> meta })                    // keep the input's meta (viruses/plasmids)
            .map { items -> [items.last(), items[0..-2]] }                      // [meta, [.tbl.gz, ...]]
    }

    // Steps above are optional: hmmer_tables / amr_gff stay channel.empty() when skipped,
    // so callers must not depend on them being emitted (use .ifEmpty).
    emit:
    hmmer_tables   = hmmer_tables                // channel: [ meta, [ tbl.gz, ... ] ]
    amr_gff        = amr_gff                     // channel: [ meta, gff ]
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}
