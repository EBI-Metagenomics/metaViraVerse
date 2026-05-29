process SEPARATE_SEQUENCES {

    label 'process_low'
    tag "${meta.id}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.75':
        'quay.io/biocontainers/biopython:1.75' }"


    input:
    tuple val(meta), path(fasta)
    val pattern

    output:
    tuple val(meta), path("${meta.id}_${pattern}.fa"), emit: chosen_sequences
    path "versions.yml",                               emit: versions

    script:
    """
    separate_sequences.py \\
       --input ${fasta} \\
       --output ${meta.id}_${pattern}.fa \\
       --pattern ${pattern} \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        biopython: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('biopython').version)")
    END_VERSIONS
    """
}
