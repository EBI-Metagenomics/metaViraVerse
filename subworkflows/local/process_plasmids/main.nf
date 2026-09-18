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
include { AMR_ANNOTATION                          } from '../../ebi-metagenomics/amr_annotation'
include { ANNOTATE_PLASMID_GFF                    } from '../../../modules/local/annotate_plasmid_gff'


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

    MOBSUITE_TYPER(
        EXTRACT_CLUSTER_FILES.out.reps_fna_uncompressed,
        params.mobsuite_db,
        [], [], [], [], [], [], [],
        true,
        true
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
        EXTRACT_CLUSTER_FILES.out.reps_gff
            .join(PLASQUID_WORKFLOW.out.protein_report)
            .join(MOBSUITE_TYPER.out.biomarker_report, remainder: true)
            .join(AMR_ANNOTATION.out.gff, remainder: true)
            .map { meta, gff, protein_report, biomarker_report, amr_gff ->
                tuple(meta, gff, protein_report, biomarker_report ?: [], amr_gff ?: [])
            }
    )
    ch_versions = ch_versions.mix(ANNOTATE_PLASMID_GFF.out.versions)

    emit:

    clustering_tsv        = CLUSTERING.out.clusters_tsv  // [meta, tsv]
    reps_tsv              = EXTRACT_CLUSTER_FILES.out.reps_list
    reps_seqs             = EXTRACT_CLUSTER_FILES.out.reps_fna_compressed      // compressed
    reps_proteins         = EXTRACT_CLUSTER_FILES.out.reps_faa_compressed      // compressed
    reps_gff              = EXTRACT_CLUSTER_FILES.out.reps_gff
    reps_gff_plasquid     = ANNOTATE_PLASMID_GFF.out.gff                        // reps_gff + plaSquid/MOB-suite/AMR evidence attributes
    versions              = ch_versions                 // channel: [ path(versions.yml) ]

}
