process GENERATE_TAXONOMY_TABLE {

    label 'process_low'
    tag "${meta.id}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.75':
        'quay.io/biocontainers/biopython:1.75' }"

    input:
      tuple val(meta), path(table)
      tuple val(meta2), path(metadata_file)

    output:
      tuple val(meta), path("*_taxonomy_counts.tsv"), path("*_metadata_table.tsv"), emit: taxonomy_and_metadata
      tuple val("${task.process}"), val('python'), eval("python --version 2>&1 | sed 's/Python //g'"), emit: versions_python, topic: versions

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
}
