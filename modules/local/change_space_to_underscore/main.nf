nextflow.enable.types = true

process CHANGE_SPACE_TO_UNDERSCORE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' &&!task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bash:5.2.37--06dbc4169cb39ae0' :
        'community.wave.seqera.io/library/bash:5.2.37--ae00789afb795adf' }"

    input:
    record(
        meta: Map,
        fasta: Path
    )

    output:
    record(
        meta: meta,
        sanitised_fasta: file("*_sanitised.fasta")
    )

    topic:
    tuple(task.process, 'bash', eval("bash --version | head -n1 | sed 's/.*version //; s/ .*//'")) >> 'versions'

    script:
    """
    sed 's/ /_/g' ${fasta} > ${fasta.baseName}_sanitised.fasta
    """

    stub:
    """
    touch ${meta.id}_sanitised.fasta
    """
}
