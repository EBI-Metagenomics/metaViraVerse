/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { CLUSTERING                              } from '../clustering/main'
include { EXTRACT_CLUSTER_FILES                   } from '../extract_cluster_files/main'
include { INDEX_RESULTS                           } from '../index_results/main'
include { PLASQUID_WORKFLOW                       } from '../plasquid_workflow/main'
include { MOBSUITE_TYPER                          } from '../../../modules/nf-core/mobsuite/typer/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PROCESS_PLASMIDS {

    take:
    sequences
    combined_faa
    combined_gff
    mapfile

    main:

    ch_versions = channel.empty()

    //
    // Cluster sequences
    //
    sequences
        .filter { _meta, seqs ->
            seqs && seqs.size() > 0
        }
        .set { samples_seqs }

    CLUSTERING(
       samples_seqs,
       'gani',
       false,
       0.35,
       false
    )
    ch_versions = ch_versions.mix(CLUSTERING.out.versions)

    //
    // ----------- Extract cluster files: FNA, FAA, TSV, GFF -----------
    //
    EXTRACT_CLUSTER_FILES (
       CLUSTERING.out.clusters_tsv,
       combined_gff,
       mapfile,
       sequences,
       combined_faa,
       'plasmids'
    )
    ch_versions = ch_versions.mix(EXTRACT_CLUSTER_FILES.out.versions)

    //
    // ----------- post-processing (indexing) for website -----------
    //

    INDEX_RESULTS (
        EXTRACT_CLUSTER_FILES.out.reps_gff,
        EXTRACT_CLUSTER_FILES.out.reps_fna_uncompressed,
        EXTRACT_CLUSTER_FILES.out.reps_faa_uncompressed,
        false
    )

    //
    // ----------- Analysis -----------
    //

    PLASQUID_WORKFLOW(
        EXTRACT_CLUSTER_FILES.out.reps_fna_uncompressed,
        EXTRACT_CLUSTER_FILES.out.reps_faa_uncompressed
    )

    MOBSUITE_TYPER(
        EXTRACT_CLUSTER_FILES.out.reps_fna_uncompressed,
        params.mobsuite_db,
        [], [], [], [], [], [], [],
        true,
        true
    )

    emit:

    clustering_tsv = CLUSTERING.out.clusters_tsv  // [meta, tsv]
    reps_tsv       = EXTRACT_CLUSTER_FILES.out.reps_list
    reps_seqs      = EXTRACT_CLUSTER_FILES.out.reps_fna_compressed      // compressed
    reps_proteins  = EXTRACT_CLUSTER_FILES.out.reps_faa_compressed      // compressed
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}
