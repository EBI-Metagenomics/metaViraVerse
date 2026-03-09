process BACPHLIP {

    label 'process_low'
    tag "combined"
    container "quay.io/microbiome-informatics/bacphlip:v0.9.3-alpha"

    input:
    tuple val(meta), path(fasta)

    output:
    path("*.bacphlip"),        emit: bacphlip_table
    path "versions.yml",       emit: versions

    script:
    """
    bacphlip \\
        --multi_fasta \\
        -f \\
        --input_file ${fasta}

    mv "${fasta}.bacphlip" ${meta.id}.bacphlip

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bacphlip: 0.9.3-alpha

    END_VERSIONS
    """
}
