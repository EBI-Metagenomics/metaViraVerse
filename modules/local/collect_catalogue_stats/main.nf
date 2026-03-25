/*
 * module to generate JSON for main website page
*/
process COLLECT_CATALOGUE_STATS {

    label 'process_low'
    tag "${meta.id}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.75':
        'quay.io/biocontainers/biopython:1.75' }"

    input:
    tuple val(meta), path(viral_seqs)
    tuple val(meta), path(plasmids)
    tuple val(meta), path(prophages)
    path(metadata)
    path(clusters_viruses)
    path(clusters_plasmids)

    output:
    tuple val(meta), path("*.json"),     emit: catalogue_json
    path "versions.yml",                 emit: versions

    script:
    """
    collect_catalogue_stats.py \\
      --viral-sequences ${viral_seqs} \\
      --plasmids ${plasmids} \\
      --prophages ${prophages} \\
      --metadata ${metadata} \\
      --clusters-viruses ${clusters_viruses} \\
      --clusters-plasmids ${clusters_plasmids} \\
      --proteins-viruses
      --proteins-plasmids

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        biopython: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('biopython').version)")
    END_VERSIONS
    """
}
