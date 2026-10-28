nextflow.enable.types = true

process PLOT_ITOL {

    label 'process_low'
    tag "${meta.id}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.75':
        'quay.io/biocontainers/biopython:1.75' }"

    input:
    record(
        meta: Map,
        table: Path,
        metadata: Path
    )

    output:
    record(
        meta: meta,
        tree_nwk: file("*_tree.nwk"),
        tree_labels: file("*_labels.txt", optional: true),
        tree_counts: file("*_counts.txt"),
        tree_types: file("*.type.binary.txt", optional: true),     // genome/assembly
        tree_biomes: file("*.biomes.binary.txt", optional: true)
    )

    topic:
    tuple(task.process, 'python', eval("python --version 2>&1 | sed 's/Python //g'")) >> 'versions'

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
    """

    stub:
    """
    echo stub > stub_tree.nwk
    echo stub > stub_labels.txt
    echo stub > stub_counts.txt
    echo stub > stub.type.binary.txt
    echo stub > stub.biomes.binary.txt
    """
}
