nextflow.enable.types = true

process GENERATE_TAXONOMY_TABLE {

    label 'process_low'
    tag "${meta.id}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.75':
        'quay.io/biocontainers/biopython:1.75' }"

    input:
    record(
        meta: Map,
        table: Path
    )
    metadata_file: Path

    output:
    record(
        meta: meta,
        taxonomy_counts: file("*_taxonomy_counts.tsv"),
        metadata_table: file("*_metadata_table.tsv")
    )

    topic:
    tuple(task.process, 'python', eval("python --version 2>&1 | sed 's/Python //g'")) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    generate_taxonomy_table.py \\
        ${args} \\
        --input ${table} \\
        --output ${prefix}_taxonomy_counts.tsv \\
        --meta ${metadata_file} \\
        --meta-output ${prefix}_metadata_table.tsv
    """

    stub:
    """
    echo stub > stub_taxonomy_counts.tsv
    echo stub > stub_metadata_table.tsv
    """
}
