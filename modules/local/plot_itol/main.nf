process PLOT_ITOL {

    label 'process_low'
    tag "${meta.id}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.75':
        'quay.io/biocontainers/biopython:1.75' }"

    input:
      tuple val(meta), path(table), path(metadata)

    output:
      tuple val(meta), path("*_tree.nwk"),   emit: tree_nwk
      tuple val(meta), path("*_labels.txt"), emit: tree_labels
      tuple val(meta), path("*_counts.txt"), emit: tree_counts
      path "versions.yml",                   emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    plot_itol.py \\
        ${args} \\
        --input ${table} \\
        --meta ${metadata} \\
        --output ${prefix}_tree.nwk \\
        --labels ${prefix}_labels.txt \\
        --counts ${prefix}_counts.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
