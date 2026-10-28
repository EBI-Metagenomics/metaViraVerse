nextflow.enable.types = true

process SUMMARISE_ANNOTATIONS {

    label 'process_low'
    tag "${meta.id}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.75':
        'quay.io/biocontainers/biopython:1.75' }"


    input:
    record(
        meta: Map,
        tbl: Path,
        metadata: Path
    )

    output:
    record(
        meta: meta,
        summary_tsv: file("${meta.id}_summary.tsv.gz")
    )

    topic:
    tuple(task.process, 'python', eval("python --version 2>&1 | sed 's/Python //g'")) >> 'versions'

    script:
    """
    summarise_annotations.py \\
       -i ${tbl} \\
       -m ${metadata} \\
       --output ${meta.id}_summary.tsv \\
       --compress
    """
}
