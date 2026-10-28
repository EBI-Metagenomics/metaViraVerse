nextflow.enable.types = true

process ANICALC {

    label 'process_low'
    tag "$meta.id"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/checkv:1.0.3--pyhdfd78af_0':
        'biocontainers/checkv:1.0.3--pyhdfd78af_0' }"

    input:
    record(
        meta: Map,
        tsv: Path
    )

    output:
    record(
        meta: meta,
        calculated_ani_tsv: file("${meta.id}_pairani.tsv")
    )

    topic:
    tuple(task.process, 'checkv', eval("checkv -h 2>&1  | sed -n 's/^.*CheckV v//; s/: assessing.*//; 1p'")) >> 'versions'

    script:
    """
    anicalc.py \\
       -i ${tsv} \\
       -o ${meta.id}_pairani.tsv
    """
}
