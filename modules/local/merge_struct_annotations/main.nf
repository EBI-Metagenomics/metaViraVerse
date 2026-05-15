/*
 * modules/local/merge_struct_annotations/main.nf
 *
 * Merges all structural annotation TSVs into the existing
 * viral_sequences_reps_stats.tsv, appending 10 new structural columns.
 *
 * Input files:
 *   reps_stats      — existing pipeline stats TSV (first column = seq_id)
 *   confidence_tsv  — ESMFold per-sequence pLDDT + status
 *   bfvd_hits       — Foldseek BFVD best-hit TSV
 *   ecod_tsv        — ECOD domain annotation TSV (with status column)
 *
 * New columns appended:
 *   struct_predicted, mean_plddt, struct_status,
 *   bfvd_top_hit, bfvd_tm_score, bfvd_evalue, bfvd_pident,
 *   ecod_xgroup, ecod_hgroup, ecod_fgroup
 *
 * Design: every input sequence gets a row even if prediction failed
 * or pLDDT was too low — status column provides the audit trail.
 */

process MERGE_STRUCT_ANNOT {
    tag "${meta.id}"
    label 'process_low'

    conda 'conda-forge::python=3.10'
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/python:3.10' :
        'python:3.10-slim' }"

    input:
    tuple val(meta), path(reps_stats)
    tuple val(meta), path(confidence_tsv)
    tuple val(meta), path(bfvd_hits)
    tuple val(meta), path(ecod_tsv)

    output:
    tuple val(meta), path("${meta.id}_struct_stats.tsv"), emit: struct_stats_tsv

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python3 ${projectDir}/scripts/merge_struct_annotations.py \\
        --reps-stats "${reps_stats}"      \\
        --confidence  "${confidence_tsv}" \\
        --bfvd-hits   "${bfvd_hits}"      \\
        --ecod        "${ecod_tsv}"       \\
        --output      "${meta.id}_struct_stats.tsv"
    """

    stub:
    """
    # Copy header from reps_stats and add structural columns
    head -1 "${reps_stats}" | tr -d '\n' > "${meta.id}_struct_stats.tsv"
    printf "\tstruct_predicted\tmean_plddt\tstruct_status\t"  >> "${meta.id}_struct_stats.tsv"
    printf "bfvd_top_hit\tbfvd_tm_score\tbfvd_evalue\t"      >> "${meta.id}_struct_stats.tsv"
    printf "bfvd_pident\tecod_xgroup\tecod_hgroup\tecod_fgroup\n" >> "${meta.id}_struct_stats.tsv"
    """
}
