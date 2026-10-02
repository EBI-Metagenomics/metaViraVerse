/*
 * module to generate JSON for main website page
*/
process COLLECT_CATALOGUE_STATS {

    label 'process_low'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.75':
        'quay.io/biocontainers/biopython:1.75' }"

    input:
    tuple val(meta1), path(json)
    tuple val(meta4), path(metadata)
    tuple val(meta5), path(clusters_viruses)
    tuple val(meta6), path(clusters_plasmids)
    tuple val(meta9), path(excluded_metadata)
    tuple val(meta10), path(initial_metadata)

    output:
    path("*.json"),                      emit: catalogue_json
    path "versions.yml",                 emit: versions

    script:
    """
    collect_catalogue_stats.py \\
      --initial-json ${json} \\
      --metadata ${metadata} \\
      --clusters-viruses ${clusters_viruses} \\
      --clusters-plasmids ${clusters_plasmids} \\
      --excluded-metadata ${excluded_metadata} \\
      --initial-metadata ${initial_metadata} \\
      -o catalogue.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        biopython: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('biopython').version)")
    END_VERSIONS
    """
}
