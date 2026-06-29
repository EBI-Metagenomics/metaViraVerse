/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FIND_CONCATENATE as CONCATENATE_BACPHLIP     } from '../../modules/nf-core/find/concatenate'
include { SEQKIT_SPLIT2 as CHUNK_FNA                   } from '../../modules/nf-core/seqkit/split2'

include { BACPHLIP                                     } from '../../modules/local/bacphlip'
include { BUILD_FINAL_GFF                              } from '../../modules/local/build_final_gff'

include { EXTRACT_CLUSTER_FILES                        } from './extract_cluster_files'
include { CLUSTERING                                   } from './clustering'
include { HOST_DETECTION                               } from './host_detection_subwf'
include { INDEX_RESULTS                                } from './index_results'
include { PROTEINS_PROCESSING                          } from './proteins_subwf'
include { TAXONOMY_ASSIGNMENT                          } from './taxonomy_subwf'

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

    ch_versions = channel.empty()

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
    ch_versions = ch_versions.mix(CLUSTERING.out.versions)

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
    ch_versions = ch_versions.mix(EXTRACT_CLUSTER_FILES.out.versions)

    //
    // ----------- Chunk nucleotide sequences for cluster reps -----------
    //

    CHUNK_FNA (
        EXTRACT_CLUSTER_FILES.out.reps_fna_uncompressed,
        [],                                        // length: (disabled) max number of nucleotides per chunk
        params.nucleotide_fasta_chunksize,         // size: max number of sequences per chunk
    )
    ch_versions = ch_versions.mix(CHUNK_FNA.out.versions)
    def ch_fna_chunks = CHUNK_FNA.out.chunked_output.transpose()

    //
    // ----------- Host assignment -----------
    // IPHOP process has chunking with another fasta size
    //

    HOST_DETECTION (
        EXTRACT_CLUSTER_FILES.out.reps_fna_uncompressed,
        params.catalogues_metadata
    )
    ch_versions = ch_versions.mix(HOST_DETECTION.out.versions)

    //
    // ----------- Lifestyle -----------
    // predicting bacteriophage lifestyle from conserved protein domains
    // running on chunked FNA file
    //

    BACPHLIP(
        ch_fna_chunks
    )
    ch_versions = ch_versions.mix(BACPHLIP.out.versions)

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
       params.additional_model_data,
       params.ncbi_db,
       params.factor_file ? params.factor_file: channel.value(false)
    )
    ch_versions = ch_versions.mix(TAXONOMY_ASSIGNMENT.out.versions)

    //
    // ----------- Proteins processing -----------
    //

    PROTEINS_PROCESSING(
       EXTRACT_CLUSTER_FILES.out.reps_faa_uncompressed.join(EXTRACT_CLUSTER_FILES.out.reps_gff)
    )
    ch_versions = ch_versions.mix(PROTEINS_PROCESSING.out.versions)

    //
    // ----------- Build final aggregated GFF -----------
    //

    BUILD_FINAL_GFF (
        EXTRACT_CLUSTER_FILES.out.reps_gff,
        CONCATENATE_BACPHLIP.out.file_out,
        PROTEINS_PROCESSING.out.hmmer_tables,
        PROTEINS_PROCESSING.out.amr_gff
    )
    ch_versions = ch_versions.mix(BUILD_FINAL_GFF.out.versions)

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

    clustering_tsv = CLUSTERING.out.clusters_tsv  // [meta, tsv]
    reps_tsv       = EXTRACT_CLUSTER_FILES.out.reps_list
    reps_seqs      = EXTRACT_CLUSTER_FILES.out.reps_fna_compressed  // compressed
    reps_proteins  = EXTRACT_CLUSTER_FILES.out.reps_faa_compressed  // compressed
    vitap_best     = TAXONOMY_ASSIGNMENT.out.vitap_best
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}
