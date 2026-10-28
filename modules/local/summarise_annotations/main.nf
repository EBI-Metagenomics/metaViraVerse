process SUMMARISE_ANNOTATIONS {

    label 'process_low'
    tag "${meta.id}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.75':
        'quay.io/biocontainers/biopython:1.75' }"


    input:
    tuple val(meta), path(tbl), path(metadata)

    output:
    tuple val(meta), path("${meta.id}_summary.tsv.gz"), emit: summary_tsv
    tuple val("${task.process}"), val('python'), eval("python --version 2>&1 | sed 's/Python //g'"), emit: versions_python, topic: versions

    script:
    """
    summarise_annotations.py \\
       -i ${tbl} \\
       -m ${metadata} \\
       --output ${meta.id}_summary.tsv \\
       --compress
    """
}
