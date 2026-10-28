nextflow.enable.types = true

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

// params of this workflow and of the subworkflows it calls
record PlasmidParams {
    cluster_vclust: Boolean
    nucleotide_fasta_chunksize: Integer
    mobsuite_db: Path?
    skip_mmseqs: Boolean
    skip_amrfinderplus_for_plasmids: Boolean
    skip_deeparg_for_plasmids: Boolean
    skip_rgi_for_plasmids: Boolean
    // plaSquid (see PlasquidParams)
    repsearch_db: Path?
    repfilter_db: Path?
    rna_inc_db: Path?
    incsearch_db: Path?
    mobsearch_db: Path?
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

record PlasmidSequences {
    meta: Map
    fna: Path
}

record MobsuiteReports {
    meta: Map
    biomarker_report: Path?   // null when MOB-suite found no biomarkers
    mob_report: Path?
}

record PlasmidCatalogue {
    meta: Map
    clusters_tsv: Path
    reps_list: Path
    fna_gz: Path               // representative sequences, compressed
    faa_gz: Path               // representative proteins, compressed
    reps_gff: Path
    annotated_gff: Path        // reps_gff + plaSquid/MOB-suite/AMR evidence attributes
    mobility_stats: Path       // plaSquid record counts + MOB-suite conjugative/mobilizable/non_mobilizable counts
}

workflow PROCESS_PLASMIDS {

    take:
    sequences: Channel<PlasmidSequences>
    combined_faa: Value<Path>
    combined_gff: Value<Path>
    opts: PlasmidParams

    main:

    //
    // --------- Cluster DNA sequences
    //
    samples_seqs = sequences.filter { r -> r.fna && r.fna.size() > 0 }

    clusters = CLUSTERING(
       samples_seqs,
       'gani',
       0,
       0.35,
       0,
       opts.cluster_vclust
    )

    //
    // ----------- Extract cluster files: FNA, FAA, TSV, GFF -----------
    //
    reps = EXTRACT_CLUSTER_FILES (
       clusters,
       sequences,
       combined_gff,
       combined_faa,
       null,
       'plasmids'
    )

    //
    // ----------- post-processing (indexing) for website -----------
    //

    INDEX_RESULTS (
        reps.map { r -> record(meta: r.meta, reps_gff: r.reps_gff, fna: r.fna, faa: r.faa) },
        channel.empty()
    )

    //
    // ----------- plasquid analysis -----------
    //

    plasquid = PLASQUID_WORKFLOW(
        reps.map { r -> record(meta: r.meta, fna: r.fna, faa: r.faa, gff: r.reps_gff) },
        opts
    )


    //
    // ----------- MOB-suite typer run on chunked fasta
    //

    chunks = CHUNK_FNA (
        reps.map { r -> tuple(r.meta, r.fna) },
        [],                                        // length: (disabled) max number of nucleotides per chunk
        opts.nucleotide_fasta_chunksize,     // size: max number of sequences per chunk
    )
    // a single chunk is emitted as a path, several as a list: transpose handles both
    def ch_fna_chunks = chunks.chunked_output.transpose()

    mobsuite = MOBSUITE_TYPER(
        ch_fna_chunks,
        opts.mobsuite_db,
        [], [], [], [], [], [], [],
        true,
        true
    )

    biomarker_report = CONCATENATE_MOBSUITE_BIOMARKER_REPORT (
        mobsuite.biomarker_report.groupTuple(),
        1
    )

    CONCATENATE_MOBSUITE_MGE_REPORT (
        mobsuite.mge_report.groupTuple(),
        1
    )

    mob_report = CONCATENATE_MOBSUITE_REPORT (
        mobsuite.report.groupTuple(),
        1
    )

    // MOB-suite reports are `optional: true`: join with remainder so a missing
    // file doesn't stall the processes below; its field is then null
    mobsuite_reports = reps
        .map { r -> record(meta: r.meta) }
        .join(biomarker_report.file_out.map { meta, f -> record(meta: meta, biomarker_report: f) }, by: 'meta', remainder: true)
        .join(mob_report.file_out.map { meta, f -> record(meta: meta, mob_report: f) }, by: 'meta', remainder: true)
        .map { r -> record(meta: r.meta, biomarker_report: r.biomarker_report, mob_report: r.mob_report) }


    //
    // ----------- Combine plaSquid and MOB-suite mobility counts -----------
    //
    // plaSquid's protein_report record counts + MOB-suite's predicted_mobility counts.

    mobility = UPDATE_MOBILITY_STATS(
        plasquid.join(mobsuite_reports, by: 'meta')
            .map { r -> record(meta: r.meta, mobility_json: r.mobility_stats, mob_report: r.mob_report) }
    )


    //
    // ----------- Proteins processing -----------
    // AMR detection, hmmer annotation - skipped, protein clustering with mmseqs
    //

    proteins = PROTEINS_PROCESSING(
       reps.map { r -> record(meta: r.meta, faa: r.faa, gff: r.reps_gff) },
       opts.skip_amrfinderplus_for_plasmids,
       opts.skip_deeparg_for_plasmids,
       opts.skip_rgi_for_plasmids,
       true,
       opts.skip_mmseqs,
       null,
       opts
    )


    //
    // ----------- Add plaSquid RIP/MOB/Inc, MOB-suite biomarker and AMR evidence to the representative GFF -----------
    //
    // Missing evidence (null field) is treated by the script as "no evidence supplied".

    annotated_gff = ANNOTATE_PLASMID_GFF(
        reps.join(plasquid, by: 'meta')
            .join(mobsuite_reports, by: 'meta')
            .join(proteins, by: 'meta')
            .map { r ->
                record(
                    meta: r.meta,
                    gff: r.reps_gff,
                    protein_report: r.protein_report,
                    biomarker_report: r.biomarker_report,
                    mob_report: r.mob_report,
                    amr_gff: r.amr_gff
                )
            }
    )

    catalogue = clusters
        .join(reps, by: 'meta')
        .join(annotated_gff.map { r -> record(meta: r.meta, annotated_gff: r.gff) }, by: 'meta')
        .join(mobility, by: 'meta')
        .map { r ->
            record(
                meta: r.meta, clusters_tsv: r.clusters_tsv, reps_list: r.reps_list,
                fna_gz: r.fna_gz, faa_gz: r.faa_gz, reps_gff: r.reps_gff,
                annotated_gff: r.annotated_gff, mobility_stats: r.mobility_stats
            )
        }

    emit:
    catalogue: Channel<PlasmidCatalogue> = catalogue
}
