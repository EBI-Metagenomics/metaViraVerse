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

workflow PROCESS_VIRAL_SEQUENCES {

    take:
    sequences
    gff
    faa
    combined_metadata
    mapfile

    main:

    //
    // ----------- Cluster sequences -----------
    //

    sequences
        .filter { meta, seqs ->
            seqs && seqs.size() > 0
        }
        .set { samples_seqs }

    CLUSTERING(
       samples_seqs,
       'ani',
       false,
       false,
       0.95
    )

    //
    // ----------- Extract cluster files: FNA, FAA, TSV, GFF -----------
    //
    EXTRACT_CLUSTER_FILES (
       CLUSTERING.out.clusters_tsv,
       gff,
       mapfile,
       sequences,
       faa,
       'viruses'
    )

    //
    // ----------- Chunk nucleotide sequences for cluster reps -----------
    //

    CHUNK_FNA (
        EXTRACT_CLUSTER_FILES.out.reps_fna_uncompressed,
        [],                                        // length: (disabled) max number of nucleotides per chunk
        params.nucleotide_fasta_chunksize,         // size: max number of sequences per chunk
    )
    def ch_fna_chunks = CHUNK_FNA.out.chunked_output.transpose()

    //
    // ----------- Host assignment -----------
    // IPHOP process has chunking with another fasta size
    //

    HOST_DETECTION (
        EXTRACT_CLUSTER_FILES.out.reps_fna_uncompressed,
        params.custom_genomes_metadata,
        mapfile
    )

    //
    // ----------- Lifestyle -----------
    // predicting bacteriophage lifestyle from conserved protein domains
    // running on chunked FNA file
    //

    BACPHLIP(
        ch_fna_chunks
    )

    CONCATENATE_BACPHLIP (
        BACPHLIP.out.bacphlip_table.groupTuple(),
        1
    )

    //
    // ----------- Taxonomy ViPhOGs and VITAP -----------
    //

    TAXONOMY_ASSIGNMENT (
       combined_metadata,
       ch_fna_chunks,
       EXTRACT_CLUSTER_FILES.out.reps_faa_uncompressed,
       EXTRACT_CLUSTER_FILES.out.reps_gff,
       params.skip_vitap,
       params.skip_genomad,
       params.skip_viphogs,
       params.viphog_db,
       params.additional_viphog_model_data,
       params.ncbi_db,
       params.viphogs_factor_file ? params.viphogs_factor_file: []
    )

    //
    // ----------- Proteins processing -----------
    // AMR detection, hmmer annotation, protein clustering with phammseqs
    //

    PROTEINS_PROCESSING(
       EXTRACT_CLUSTER_FILES.out.reps_faa_uncompressed.join(EXTRACT_CLUSTER_FILES.out.reps_gff),
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
        EXTRACT_CLUSTER_FILES.out.reps_gff,
        CONCATENATE_BACPHLIP.out.file_out,
        // optional: empty when --annotation_db is not set / all AMR tools are skipped
        PROTEINS_PROCESSING.out.hmmer_tables.ifEmpty([[:], []]),
        PROTEINS_PROCESSING.out.amr_gff.ifEmpty([[:], []])
    )

    //
    // ----------- post-processing (indexing) for website -----------
    //

    INDEX_RESULTS (
        EXTRACT_CLUSTER_FILES.out.reps_gff,
        EXTRACT_CLUSTER_FILES.out.reps_fna_uncompressed,
        EXTRACT_CLUSTER_FILES.out.reps_faa_uncompressed,
        CONCATENATE_BACPHLIP.out.file_out
    )

    emit:

    clustering_tsv      = CLUSTERING.out.clusters_tsv  // [meta, tsv]
    reps_tsv            = EXTRACT_CLUSTER_FILES.out.reps_list
    reps_seqs           = EXTRACT_CLUSTER_FILES.out.reps_fna_compressed  // compressed
    reps_proteins       = EXTRACT_CLUSTER_FILES.out.reps_faa_compressed  // compressed
    vitap_best          = TAXONOMY_ASSIGNMENT.out.vitap_best
    viphogs_assign      = TAXONOMY_ASSIGNMENT.out.viphogs_taxonomy
    genomad_assign      = TAXONOMY_ASSIGNMENT.out.genomad_taxonomy
    iphop_host_genome   = HOST_DETECTION.out.iphop_host_genome
    iphop_host_genus    = HOST_DETECTION.out.iphop_host_genus
    spacepharer_host    = HOST_DETECTION.out.spacepharer_host

}
