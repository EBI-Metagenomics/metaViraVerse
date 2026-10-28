nextflow.enable.types = true

process BACPHLIP {

    label 'process_medium'
    tag "${fasta.baseName}"
    container "quay.io/microbiome-informatics/bacphlip:v0.9.3-alpha"

    input:
    record(
        meta: Map,
        fasta: Path
    )

    output:
    record(
        meta: meta,
        bacphlip_table: file("*.bacphlip")
    )

    topic:
    tuple(task.process, 'bacphlip', '0.9.3-alpha') >> 'versions'

    script:
    """
    bacphlip \\
        --multi_fasta \\
        -f \\
        --input_file ${fasta}
    """
}
