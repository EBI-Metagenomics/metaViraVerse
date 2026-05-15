/*
 * modules/local/merge_struct_annotations/main.nf
 *
 * Merges all structural annotation TSVs into the existing
 * viral_sequences_reps_stats.tsv, adding 9 new structural columns.
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
        --reps-stats "${reps_stats}" \\
        --confidence  "${confidence_tsv}" \\
        --bfvd-hits   "${bfvd_hits}" \\
        --ecod        "${ecod_tsv}" \\
        --output      "${meta.id}_struct_stats.tsv"
    """

    stub:
    """
    touch "${meta.id}_struct_stats.tsv"
    """
}
