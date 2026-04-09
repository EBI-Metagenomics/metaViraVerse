include { PHAMMSEQS                           } from '../../modules/local/phammseqs'
include { SUMMARISE_ANNOTATIONS               } from '../../modules/local/summarise_annotations'

include { HMMER_HMMSEARCH                     } from '../../modules/nf-core/hmmer/hmmsearch'
include { SEQKIT_SPLIT2                       } from '../../modules/nf-core/seqkit/split2/main'
include { CAT_CAT as CONCATENATE_HMMER_TBLOUT } from '../../modules/nf-core/cat/cat/main'

include { AMR_ANNOTATION                      } from '../ebi-metagenomics/amr_annotation'


workflow PROTEINS_PROCESSING {

    take:
    input  // (meta, reps_faa.uncompressed, reps_gff)

    main:

    ch_versions = channel.empty()

    ch_proteins = input.map{ meta, faa, _gff -> tuple(meta, faa) }

    /*********************/
    /* PROTEINS CHUNKING */
    /*********************/

    // Chunk protein FASTA into files with N sequences each for parallel annotation
    // Default: 50,000 sequences per chunk (override with --protein_annotation_fasta_chunksize)
    SEQKIT_SPLIT2(
        COMBINED_GENE_CALLER.out.faa,
        [],                                        // length: (disabled) max number of nucleotides per chunk
        params.protein_annotation_fasta_chunksize, // size: max number of sequences per chunk
    )
    ch_versions = ch_versions.mix(SEQKIT_SPLIT2.out.versions)

    def ch_protein_chunks = SEQKIT_SPLIT2.out.assembly.transpose()

    //
    // -------- Antimicrobial resistence detection
    //
    AMR_ANNOTATION (
        input,
        params.amrfinderplus_db,
        params.deeparg_db,
        params.deeparg_db_version,
        params.deeparg_model,
        params.deeparg_tool_version,
        params.rgi_db,
        params.skip_amrfinderplus,
        params.skip_deeparg,
        params.skip_rgi
    )

    if (params.phammseqs) {
        //
        // -------- Assort phage protein sequences into phamilies using MMseqs2
        //
        PHAMMSEQS(
           ch_proteins
        )
        ch_versions = ch_versions.mix(PHAMMSEQS.out.versions)
        // TODO: maybe run PHAMCLUST on small dataset
    }

    //
    // -------- Functional annotation with HMMER
    //
    def hmm_ch = channel
        .fromPath("${params.annotation_db}/*.hmm.gz")
        .map { hmm_file -> tuple(hmm_file.baseName, hmm_file) }

    ch_protein_chunks.view()
    hmm_ch.combine(ch_protein_chunks).view()

    def hmmsearch_input = hmm_ch
        .combine(ch_proteins)
        .map { hmm_id, hmm_file, faa_id, faa_file ->
            def meta = [
                id: hmm_id.replace('.hmm', '')
            ]
            tuple(meta, hmm_file, faa_file, true, true, true)
        }
    HMMER_HMMSEARCH (
        hmmsearch_input
    )
    ch_versions = ch_versions.mix(HMMER_HMMSEARCH.out.versions)

    //
    // ----- Add metadata and generate a full summary
    //
    def hmm_metadata = channel
        .fromPath("${params.annotation_db}/*.tsv")
        .map { hmm_file ->
             def meta = [
                id: hmm_file.baseName
             ]
             tuple(meta, hmm_file)
        }

    SUMMARISE_ANNOTATIONS(
        HMMER_HMMSEARCH.out.target_summary.join(hmm_metadata)
    )
    ch_versions = ch_versions.mix(SUMMARISE_ANNOTATIONS.out.versions)

    emit:
    hmmer_tables   = HMMER_HMMSEARCH.out.target_summary
                      .map { meta, table -> table }                  // Extract just the table files
                      .collect()                                     // Gather all tables into a list
                      .map { tables -> [[id: 'viruses'], tables] }   // [meta, [.tbl.gz, ...]]
    amr_gff        = AMR_ANNOTATION.out.gff
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}