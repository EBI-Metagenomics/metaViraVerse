nextflow.enable.types = true

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FIND_CONCATENATE as CONCATENATE_BACPHLIP     } from '../../../modules/nf-core/find/concatenate'
include { SEQKIT_SPLIT2 as CHUNK_FNA                   } from '../../../modules/nf-core/seqkit/split2'

include { BACPHLIP                                     } from '../../../modules/local/bacphlip'
include { BUILD_FINAL_GFF                              } from '../../../modules/local/build_final_gff'

include { EXTRACT_CLUSTER_FILES                        } from '../extract_cluster_files/main'
include { CLUSTERING                                   } from '../clustering/main'
include { HOST_DETECTION ; HostResult                  } from '../host_detection_subwf/main'
include { INDEX_RESULTS ; IndexedReps                  } from '../index_results/main'
include { PROTEINS_PROCESSING ; ProteinsResult         } from '../proteins_subwf/main'
include { TAXONOMY_ASSIGNMENT ; TaxonomyResult         } from '../taxonomy_subwf/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// params of this workflow and of the subworkflows it calls
record ViralParams {
    cluster_vclust: Boolean
    nucleotide_fasta_chunksize: Integer
    skip_phammseqs: Boolean
    skip_amrfinderplus_for_viruses: Boolean
    skip_deeparg_for_viruses: Boolean
    skip_rgi_for_viruses: Boolean
    annotation_db: Path?
    // taxonomy (see TaxonomyParams)
    skip_vitap: Boolean
    skip_genomad: Boolean
    skip_viphogs: Boolean
    vitap_db: Path?
    genomad_db: Path?
    viphog_db: Path?
    additional_viphog_model_data: Path?
    ncbi_db: Path?
    viphogs_factor_file: Path?
    viphogs_evalue: Float
    // host detection (see HostParams)
    skip_iphop: Boolean
    iphop_db: Path?
    nucleotide_fasta_chunksize_iphop: Integer
    predict_host_from_custom_spacers: Boolean
    custom_spacers_fasta: Path?
    custom_spacers_metadata: Path?
    custom_genomes_metadata: Path?
    // proteins (see ProteinsParams)
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

record ViralSequences {
    meta: Map
    fna: Path
}

record LifestyleTable {
    meta: Map
    bacphlip: Path
}

record ViralCatalogue {
    meta: Map
    clusters_tsv: Path
    ani_tsv: Path?             // CheckV ANI clustering only
    ani_clusters_tsv: Path?
    reps_list: Path
    fna_gz: Path               // representative sequences, compressed
    faa_gz: Path               // representative proteins, compressed
    reps_gff: Path
    final_gff: Path            // reps GFF with lifestyle, HMMER and AMR annotations
    index: IndexedReps         // bgzipped and indexed reps FNA/FAA/GFF and lifestyle table
    taxonomy: TaxonomyResult
    hosts: HostResult
    proteins: ProteinsResult
}

workflow PROCESS_VIRAL_SEQUENCES {

    take:
    sequences: Channel<ViralSequences>
    gff: Value<Path>
    faa: Value<Path>
    combined_metadata: Value<Path>
    mapfile: Value<Path>
    opts: ViralParams

    main:

    //
    // ----------- Cluster sequences -----------
    //

    samples_seqs = sequences.filter { r -> r.fna && r.fna.size() > 0 }

    clusters = CLUSTERING(
       samples_seqs,
       'ani',
       0,
       0,
       0.95,
       opts.cluster_vclust
    )

    //
    // ----------- Extract cluster files: FNA, FAA, TSV, GFF -----------
    //
    reps = EXTRACT_CLUSTER_FILES (
       clusters,
       sequences,
       gff,
       faa,
       mapfile,
       'viruses'
    )

    //
    // ----------- Chunk nucleotide sequences for cluster reps -----------
    //

    chunks = CHUNK_FNA (
        reps.map { r -> tuple(r.meta, r.fna) },
        [],                                        // length: (disabled) max number of nucleotides per chunk
        opts.nucleotide_fasta_chunksize,           // size: max number of sequences per chunk
    )
    // a single chunk is emitted as a path, several as a list: transpose handles both
    def fna_chunks = chunks.chunked_output
        .transpose()
        .map { meta, fasta -> record(meta: meta, fasta: fasta) }

    //
    // ----------- Host assignment -----------
    // IPHOP process has chunking with another fasta size
    //

    hosts = HOST_DETECTION (
        reps.map { r -> record(meta: r.meta, fna: r.fna) },
        mapfile,
        opts
    )

    //
    // ----------- Lifestyle -----------
    // predicting bacteriophage lifestyle from conserved protein domains
    // running on chunked FNA file
    //

    bacphlip = BACPHLIP(
        fna_chunks
    )

    bacphlip_concat = CONCATENATE_BACPHLIP (
        bacphlip.map { r -> tuple(r.meta, r.bacphlip_table) }.groupTuple(),
        1
    )
    def bacphlip_table: Channel<LifestyleTable> = bacphlip_concat.file_out.map { meta, table -> record(meta: meta, bacphlip: table) }

    //
    // ----------- Taxonomy ViPhOGs and VITAP -----------
    //

    taxonomy = TAXONOMY_ASSIGNMENT (
       reps,
       fna_chunks,
       combined_metadata,
       opts
    )

    //
    // ----------- Proteins processing -----------
    // AMR detection, hmmer annotation, protein clustering with phammseqs
    //

    proteins = PROTEINS_PROCESSING(
       reps.map { r -> record(meta: r.meta, faa: r.faa, gff: r.reps_gff) },
       opts.skip_amrfinderplus_for_viruses,
       opts.skip_deeparg_for_viruses,
       opts.skip_rgi_for_viruses,
       opts.skip_phammseqs,
       true,
       opts.annotation_db,
       opts
    )

    //
    // ----------- Build final aggregated GFF -----------
    //

    final_gff = BUILD_FINAL_GFF (
        reps.join(bacphlip_table, by: 'meta')
            .join(proteins, by: 'meta')
            .map { r ->
                record(
                    meta: r.meta,
                    input_gff: r.reps_gff,
                    bacphlip: r.bacphlip,
                    hmmfiles: r.hmmer_tables ?: [],   // null when --annotation_db is not set
                    additional_gff: r.amr_gff         // null when all AMR tools are skipped
                )
            }
    )

    //
    // ----------- post-processing (indexing) for website -----------
    //

    indexed = INDEX_RESULTS (
        reps.map { r -> record(meta: r.meta, reps_gff: r.reps_gff, fna: r.fna, faa: r.faa) },
        bacphlip_table
    )

    catalogue = clusters
        .join(reps, by: 'meta')
        .join(final_gff, by: 'meta')
        .join(indexed.map { r -> record(meta: r.meta, index: r) }, by: 'meta')
        .join(taxonomy.map { r -> record(meta: r.meta, taxonomy: r) }, by: 'meta')
        .join(hosts.map { r -> record(meta: r.meta, hosts: r) }, by: 'meta')
        .join(proteins.map { r -> record(meta: r.meta, proteins: r) }, by: 'meta')
        .map { r ->
            record(
                meta: r.meta, clusters_tsv: r.clusters_tsv, ani_tsv: r.ani_tsv, ani_clusters_tsv: r.ani_clusters_tsv,
                reps_list: r.reps_list, fna_gz: r.fna_gz, faa_gz: r.faa_gz, reps_gff: r.reps_gff,
                final_gff: r.final_gff, index: r.index,
                taxonomy: r.taxonomy, hosts: r.hosts, proteins: r.proteins
            )
        }

    emit:
    catalogue: Channel<ViralCatalogue> = catalogue
}
