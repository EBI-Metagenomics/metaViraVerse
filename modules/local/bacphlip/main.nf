process BACPHLIP {

    label 'process_medium'
    tag "${fasta.baseName}"
    container "quay.io/microbiome-informatics/bacphlip:v0.9.3-alpha"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.bacphlip"),        emit: bacphlip_table
    tuple val("${task.process}"), val('bacphlip'), val('0.9.3-alpha'), emit: versions_bacphlip, topic: versions

    script:
    """
    bacphlip \\
        --multi_fasta \\
        -f \\
        --input_file ${fasta}
    """
}
