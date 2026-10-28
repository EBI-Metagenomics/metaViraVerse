/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { CLUSTERING                                                } from '../clustering/main'
include { EXTRACT_CLUSTER_FILES                                     } from '../extract_cluster_files/main'
include { INDEX_RESULTS                                             } from '../index_results/main'
include { PLASQUID_WORKFLOW                                         } from '../plasquid_workflow/main'
include { PROTEINS_PROCESSING                                       } from '../proteins_subwf/main'

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

    //
    // --------- Cluster DNA sequences
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
    // ----------- plasquid analysis -----------
    //

    PLASQUID_WORKFLOW(
        EXTRACT_CLUSTER_FILES.out.reps_fna_uncompressed,
        EXTRACT_CLUSTER_FILES.out.reps_faa_uncompressed,
        EXTRACT_CLUSTER_FILES.out.reps_gff
    )


    //
    // ----------- MOB-suite typer run on chunked fasta
    //

    CHUNK_FNA (
        EXTRACT_CLUSTER_FILES.out.reps_fna_uncompressed,
        [],                                        // length: (disabled) max number of nucleotides per chunk
        params.nucleotide_fasta_chunksize,   // size: max number of sequences per chunk
    )
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
    // ----------- Combine plaSquid and MOB-suite mobility counts -----------
    //
    // plaSquid's protein_report record counts + MOB-suite's predicted_mobility counts.
    // Same `optional: true` remainder-join safety as above.

    UPDATE_MOBILITY_STATS(
        PLASQUID_WORKFLOW.out.mobility_stats
            .join(CONCATENATE_MOBSUITE_REPORT.out.file_out, remainder: true)
            .map { meta, mobility_json, mob_report ->
                record(meta: meta, mobility_json: mobility_json, mob_report: mob_report)
            }
    )


    //
    // ----------- Proteins processing -----------
    // AMR detection, hmmer annotation - skipped, protein clustering with mmseqs
    //

    PROTEINS_PROCESSING(
       EXTRACT_CLUSTER_FILES.out.reps_faa_uncompressed.join(EXTRACT_CLUSTER_FILES.out.reps_gff),
       params.skip_amrfinderplus_for_plasmids,
       params.skip_deeparg_for_plasmids,
       params.skip_rgi_for_plasmids,
       true,
       params.skip_mmseqs,
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
        EXTRACT_CLUSTER_FILES.out.reps_gff.map { meta, gff -> record(meta: meta, gff: gff) },
        PLASQUID_WORKFLOW.out.protein_report.map { _meta, f -> f },
        CONCATENATE_MOBSUITE_BIOMARKER_REPORT.out.file_out.map { _meta, f -> f }.ifEmpty([]),
        CONCATENATE_MOBSUITE_REPORT.out.file_out.map { _meta, f -> f }.ifEmpty([]),
        PROTEINS_PROCESSING.out.amr_gff.map { _meta, f -> f }.ifEmpty([])
    )


    emit:

    clustering_tsv        = CLUSTERING.out.clusters_tsv  // [meta, tsv]
    reps_tsv              = EXTRACT_CLUSTER_FILES.out.reps_list
    reps_seqs             = EXTRACT_CLUSTER_FILES.out.reps_fna_compressed      // compressed
    reps_proteins         = EXTRACT_CLUSTER_FILES.out.reps_faa_compressed      // compressed
    reps_gff              = EXTRACT_CLUSTER_FILES.out.reps_gff
    reps_gff_plasquid     = ANNOTATE_PLASMID_GFF.out.map { r -> [r.meta, r.gff] }                        // reps_gff + plaSquid/MOB-suite/AMR evidence attributes
    mobility_stats        = UPDATE_MOBILITY_STATS.out.map { r -> [r.meta, r.mobility_stats] }            // plaSquid record counts + MOB-suite conjugative/mobilizable/non_mobilizable counts

}
