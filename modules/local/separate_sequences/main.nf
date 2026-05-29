process SEPARATE_SEQUENCES {

    label 'process_low'
    tag "${meta.id}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.75':
        'quay.io/biocontainers/biopython:1.75' }"


    input:
    tuple val(meta), path(fasta)
    val pattern
    path(map_file)

    output:
    tuple val(meta), path("${meta.id}_${pattern}.fasta"), emit: chosen_sequences
    path "versions.yml",                                  emit: versions

    script:
    def mapping = map_file ? "--map ${map_file}" : ""
    """
    separate_sequences.py \\
       --input ${fasta} \\
       --output ${meta.id}_${pattern}.fa \\
       --pattern ${pattern} \\
       ${mapping}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        biopython: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('biopython').version)")
    END_VERSIONS
    """
}
