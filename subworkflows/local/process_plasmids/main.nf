/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { CLUSTERING                                                } from '../clustering/main'
include { EXTRACT_CLUSTER_FILES                                     } from '../extract_cluster_files/main'
include { INDEX_RESULTS                                             } from '../index_results/main'
include { PLASQUID_WORKFLOW                                         } from '../plasquid_workflow/main'

include { AMR_ANNOTATION                                            } from '../../ebi-metagenomics/amr_annotation'

include { MOBSUITE_TYPER                                            } from '../../../modules/nf-core/mobsuite/typer/main'
include { SEQKIT_SPLIT2 as CHUNK_FNA                                } from '../../../modules/nf-core/seqkit/split2'
include { FIND_CONCATENATE as CONCATENATE_MOBSUITE_BIOMARKER_REPORT } from '../../../modules/nf-core/find/concatenate'
include { FIND_CONCATENATE as CONCATENATE_MOBSUITE_MGE_REPORT       } from '../../../modules/nf-core/find/concatenate'
include { FIND_CONCATENATE as CONCATENATE_MOBSUITE_REPORT           } from '../../../modules/nf-core/find/concatenate'

include { ANNOTATE_PLASMID_GFF                                      } from '../../../modules/local/annotate_plasmid_gff'
include { UPDATE_MOBILITY_STATS                                     } from '../../../modules/local/update_mobility_stats'


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
        EXTRACT_CLUSTER_FILES.out.reps_faa_uncompressed,
        EXTRACT_CLUSTER_FILES.out.reps_gff
    )
    ch_versions = ch_versions.mix(PLASQUID_WORKFLOW.out.versions)

    CHUNK_FNA (
        EXTRACT_CLUSTER_FILES.out.reps_fna_uncompressed,
        [],                                        // length: (disabled) max number of nucleotides per chunk
        params.nucleotide_fasta_chunksize,   // size: max number of sequences per chunk
    )
    ch_versions = ch_versions.mix(CHUNK_FNA.out.versions)
    def ch_fna_chunks = CHUNK_FNA.out.chunked_output.transpose()

    MOBSUITE_TYPER(
        ch_fna_chunks,
        params.mobsuite_db,
        [], [], [], [], [], [], [],
        true,
        true
    )

    CONCATENATE_MOBSUITE_BIOMARKER_REPORT (
        MOBSUITE_TYPER.out.biomarker_report.groupTuple(),
        1
    )

    CONCATENATE_MOBSUITE_MGE_REPORT (
        MOBSUITE_TYPER.out.mge_report.groupTuple(),
        1
    )

    CONCATENATE_MOBSUITE_REPORT (
        MOBSUITE_TYPER.out.report.groupTuple(),
        1
    )

    //
    // -------- Antimicrobial resistence detection
    //
    AMR_ANNOTATION (
        EXTRACT_CLUSTER_FILES.out.reps_faa_uncompressed.join(EXTRACT_CLUSTER_FILES.out.reps_gff),
        params.amrfinderplus_db,
        params.deeparg_db,
        params.deeparg_db_version,
        params.deeparg_model,
        params.deeparg_tool_version,
        params.rgi_db,
        false,
        false,
        false
    )

    //
    // ----------- Add plaSquid RIP/MOB/Inc, MOB-suite biomarker and AMR evidence to the representative GFF -----------
    //
    // MOBSUITE_TYPER.out.biomarker_report and AMR_ANNOTATION.out.gff (via
    // AMRINTEGRATOR) are both `optional: true` -- join with remainder so a
    // missing/never-emitted file doesn't stall the process, falling back to
    // `[]` (no file), which the script treats as "no evidence supplied" via
    // its own --biomarker-report/--amr-gff ? ... : "" checks.
    ANNOTATE_PLASMID_GFF(
        EXTRACT_CLUSTER_FILES.out.reps_gff,
        PLASQUID_WORKFLOW.out.protein_report.map { _meta, f -> f },
        CONCATENATE_MOBSUITE_BIOMARKER_REPORT.out.file_out.map { _meta, f -> f }.ifEmpty([]),
        CONCATENATE_MOBSUITE_REPORT.out.file_out.map { _meta, f -> f }.ifEmpty([]),
        AMR_ANNOTATION.out.gff.map { _meta, f -> f }.ifEmpty([])
    )
    ch_versions = ch_versions.mix(ANNOTATE_PLASMID_GFF.out.versions)

    //
    // ----------- Combine plaSquid and MOB-suite mobility counts -----------
    //
    // plaSquid's protein_report record counts + MOB-suite's predicted_mobility counts.
    // Same `optional: true` remainder-join safety as above.
    UPDATE_MOBILITY_STATS(
        PLASQUID_WORKFLOW.out.mobility_stats
            .join(CONCATENATE_MOBSUITE_REPORT.out.file_out, remainder: true)
            .map { meta, mobility_json, mob_report ->
                tuple(meta, mobility_json, mob_report ?: [])
            }
    )
    ch_versions = ch_versions.mix(UPDATE_MOBILITY_STATS.out.versions)

    emit:

    clustering_tsv        = CLUSTERING.out.clusters_tsv  // [meta, tsv]
    reps_tsv              = EXTRACT_CLUSTER_FILES.out.reps_list
    reps_seqs             = EXTRACT_CLUSTER_FILES.out.reps_fna_compressed      // compressed
    reps_proteins         = EXTRACT_CLUSTER_FILES.out.reps_faa_compressed      // compressed
    reps_gff              = EXTRACT_CLUSTER_FILES.out.reps_gff
    reps_gff_plasquid     = ANNOTATE_PLASMID_GFF.out.gff                        // reps_gff + plaSquid/MOB-suite/AMR evidence attributes
    mobility_stats        = UPDATE_MOBILITY_STATS.out.mobility_stats            // plaSquid record counts + MOB-suite conjugative/mobilizable/non_mobilizable counts
    versions              = ch_versions                 // channel: [ path(versions.yml) ]

}
