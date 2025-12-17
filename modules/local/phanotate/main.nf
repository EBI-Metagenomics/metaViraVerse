/*
 * PHANOTATE: a gene caller for phages
*/
process PHANOTATE {

    label 'process_low'
    tag "${meta.id}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/phanotate:1.6.7--ce90a2ea63fc2b28':
        'community.wave.seqera.io/library/phanotate:1.6.7--e8c9a636fb881d0e' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}.phanotate.faa"),    emit: proteins
    path "versions.yml",                                  emit: versions

    script:
    def args = task.ext.args   ?: ''
    """
    phanotate.py \\
      $args \\
      -o ${meta.id}.phanotate.faa \\
      ${fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phanotate: \$(phanotate.py --version 2>&1')
    END_VERSIONS
    """
}