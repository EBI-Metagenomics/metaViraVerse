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
include { HOST_DETECTION                               } from '../host_detection_subwf/main'
include { INDEX_RESULTS                                } from '../index_results/main'
include { PROTEINS_PROCESSING                          } from '../proteins_subwf/main'
include { TAXONOMY_ASSIGNMENT                          } from '../taxonomy_subwf/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

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
    reps_list: Path
    fna_gz: Path               // representative sequences, compressed
    faa_gz: Path               // representative proteins, compressed
    reps_gff: Path
    vitap_best: Path?
    viphogs: Path?
    genomad: Path?
    iphop_genome: Path?
    iphop_genus: Path?
    spacepharer: Path?
}

workflow PROCESS_VIRAL_SEQUENCES {

    take:
    sequences: Channel<ViralSequences>
    gff: Value<Path>
    faa: Value<Path>
    combined_metadata: Value<Path>
    mapfile: Value<Path>

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
       0.95
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
        params.nucleotide_fasta_chunksize,         // size: max number of sequences per chunk
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
        params.custom_genomes_metadata,
        mapfile
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
       params.skip_vitap,
       params.skip_genomad,
       params.skip_viphogs,
       params.viphog_db ? file(params.viphog_db) : null,
       params.additional_viphog_model_data ? file(params.additional_viphog_model_data) : null,
       params.ncbi_db ? file(params.ncbi_db) : null,
       params.viphogs_factor_file ? file(params.viphogs_factor_file) : null
    )

    //
    // ----------- Proteins processing -----------
    // AMR detection, hmmer annotation, protein clustering with phammseqs
    //

    proteins = PROTEINS_PROCESSING(
       reps.map { r -> record(meta: r.meta, faa: r.faa, gff: r.reps_gff) },
       params.skip_amrfinderplus_for_viruses,
       params.skip_deeparg_for_viruses,
       params.skip_rgi_for_viruses,
       params.skip_phammseqs,
       true,
       params.annotation_db
    )

    //
    // ----------- Build final aggregated GFF -----------
    //

    BUILD_FINAL_GFF (
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

    INDEX_RESULTS (
        reps.map { r -> record(meta: r.meta, reps_gff: r.reps_gff, fna: r.fna, faa: r.faa) },
        bacphlip_table
    )

    catalogue = clusters
        .join(reps, by: 'meta')
        .join(taxonomy, by: 'meta')
        .join(hosts, by: 'meta')
        .map { r ->
            record(
                meta: r.meta, clusters_tsv: r.clusters_tsv, reps_list: r.reps_list,
                fna_gz: r.fna_gz, faa_gz: r.faa_gz, reps_gff: r.reps_gff,
                vitap_best: r.vitap_best, viphogs: r.viphogs, genomad: r.genomad,
                iphop_genome: r.iphop_genome, iphop_genus: r.iphop_genus, spacepharer: r.spacepharer
            )
        }

    emit:
    catalogue: Channel<ViralCatalogue> = catalogue
}
