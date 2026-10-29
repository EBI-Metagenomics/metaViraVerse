nextflow.enable.types = true

process COLLECT_HOST_INFO {

    label 'process_low'
    tag "${meta.id}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.75':
        'quay.io/biocontainers/biopython:1.75' }"

    input:
    record(
        meta: Map,
        predictions: Path
    )
    spacers_metadata: Path
    mapping_file: Path
    metadata: Path

    output:
    record(
        meta: meta,
        tsv: file("custom_chosen_host.tsv"),
        lineage_comparison_json: file("*.json")
    )

    topic:
    tuple(task.process, 'python', eval("python --version 2>&1 | sed 's/Python //g'")) >> 'versions'

    script:
    """
    collect_crispr_host_info.py \\
        -p ${predictions} \\
        -c ${spacers_metadata} \\
        -m ${metadata} \\
        -o custom_chosen_host.tsv \\
        --rename-map ${mapping_file}
    """

    stub:
    """
    echo stub > custom_chosen_host.tsv
    echo stub > stub.json
    """
}
