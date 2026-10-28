nextflow.enable.types = true

/*
 * module to generate JSON for main website page
*/
process COLLECT_CATALOGUE_STATS {

    label 'process_low'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.75':
        'quay.io/biocontainers/biopython:1.75' }"

    input:
    json: Path
    metadata: Path
    clusters_viruses: Path
    clusters_plasmids: Path
    excluded_metadata: Path
    initial_metadata: Path
    clusters_viruses_prot: Path
    clusters_plasmids_prot: Path

    output:
    record(
        catalogue_json: file("*.json")
    )

    topic:
    tuple(task.process, 'python', eval("python --version 2>&1 | sed 's/Python //g'")) >> 'versions'
    tuple(task.process, 'biopython', eval("python -c \"import pkg_resources; print(pkg_resources.get_distribution('biopython').version)\"")) >> 'versions'

    script:
    """
    collect_catalogue_stats.py \\
      --initial-json ${json} \\
      --metadata ${metadata} \\
      --clusters-viruses ${clusters_viruses} \\
      --clusters-plasmids ${clusters_plasmids} \\
      --excluded-metadata ${excluded_metadata} \\
      --initial-metadata ${initial_metadata} \\
      --viruses-reps-proteins ${clusters_viruses_prot} \\
      --plasmids-reps-proteins ${clusters_plasmids_prot} \\
      -o catalogue.json
    """
}
