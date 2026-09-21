process UPDATE_MOBILITY_STATS {
    /*
     * Recompute plaSquid's mobility_stats.json using MOB-suite's biomarker report as
     * the mating-pair-formation (MPF) evidence plaSquid's own MOBSEARCH can't provide:
     * plaSquid's classification never places a contig in "conjugative" on its own
     * (relaxase-only detection), so every "mobilizable" contig that also has an MPF
     * hit in MOB-suite's plasmids_biomarker_report.txt is promoted here.
    */

    label 'process_single'
    tag "${meta.id}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.75':
        'quay.io/biocontainers/biopython:1.75' }"

    input:
    tuple val(meta), path(mobility_classification), path(mobility_json), path(biomarker_report)

    output:
    tuple val(meta), path("mobility_stats_final.json"), emit: mobility_stats
    path "versions.yml",                                emit: versions

    script:
    def biomarker_report_arg = biomarker_report ? "--biomarker-report ${biomarker_report}" : ""
    """
    update_mobility_stats.py \\
        --classification ${mobility_classification} \\
        --mobility-json ${mobility_json} \\
        ${biomarker_report_arg} \\
        --output mobility_stats_final.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
