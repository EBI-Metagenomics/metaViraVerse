nextflow.enable.types = true

process UPDATE_MOBILITY_STATS {
    /*
     * Combine plaSquid's mobility_stats.json (protein_report record counts) with
     * MOB-suite's predicted_mobility counts from the mob_typer report
     * (mobsuite_total, mobsuite_conjugative, mobsuite_mobilizable,
     * mobsuite_non_mobilizable) into one summary JSON.
    */

    label 'process_single'
    tag "${meta.id}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.75':
        'quay.io/biocontainers/biopython:1.75' }"

    input:
    record(
        meta: Map,
        mobility_json: Path,
        mob_report: Path?
    )

    output:
    record(
        meta: meta,
        mobility_stats: file("mobility_stats_final.json")
    )

    topic:
    tuple(task.process, 'python', eval("python --version 2>&1 | sed 's/Python //g'")) >> 'versions'

    script:
    def mob_report_arg = mob_report ? "--mob-report ${mob_report}" : ""
    """
    update_mobility_stats.py \\
        --mobility-json ${mobility_json} \\
        ${mob_report_arg} \\
        --output mobility_stats_final.json
    """

    stub:
    """
    echo stub > mobility_stats_final.json
    """
}
