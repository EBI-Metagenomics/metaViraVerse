process ANICLUST {

    label 'process_low'
    tag "$meta.id"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/checkv:1.0.3--pyhdfd78af_0':
        'biocontainers/checkv:1.0.3--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(ani_tsv), path(sequences)

    output:
    tuple val(meta), path("${meta.id}_clusters.tsv"), emit: clusters_tsv
    path "versions.yml"                             , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    aniclust.py \\
      --fna ${sequences} \\
      --ani ${ani_tsv} \\
      --out ${meta.id}_clusters.tsv \\
      ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkv: \$(checkv -h 2>&1  | sed -n 's/^.*CheckV v//; s/: assessing.*//; 1p')
    END_VERSIONS
    """
}
