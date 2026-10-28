process ANICALC {

    label 'process_low'
    tag "$meta.id"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/checkv:1.0.3--pyhdfd78af_0':
        'biocontainers/checkv:1.0.3--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(tsv)

    output:
    tuple val(meta), path("${meta.id}_pairani.tsv"), emit: calculated_ani_tsv
    tuple val("${task.process}"), val('checkv'), eval("checkv -h 2>&1  | sed -n 's/^.*CheckV v//; s/: assessing.*//; 1p'"), emit: versions_checkv, topic: versions

    script:
    """
    anicalc.py \\
       -i ${tsv} \\
       -o ${meta.id}_pairani.tsv
    """
}
