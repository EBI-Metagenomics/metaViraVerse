/*
 * module to generate JSON for main website page
*/
process COLLECT_CATALOGUE_STATS {

    label 'process_low'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.75':
        'quay.io/biocontainers/biopython:1.75' }"

    input:
    tuple val(meta1), path(viral_seqs)
    tuple val(meta2), path(prophages)
    tuple val(meta3), path(plasmids)
    tuple val(meta4), path(metadata)
    tuple val(meta5), path(clusters_viruses)
    path(clusters_plasmids)
    tuple val(meta6), path(viral_proteins)
    tuple val(meta7), path(plasmid_proteins)

    output:
    path("*.json"),                      emit: catalogue_json
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
      --proteins-viruses ${viral_proteins} \\
      --proteins-plasmids ${plasmid_proteins} \\
      -o catalogue.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        biopython: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('biopython').version)")
    END_VERSIONS
    """
}
