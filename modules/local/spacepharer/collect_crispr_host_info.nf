process COLLECT_HOST_INFO {

    label 'process_low'
    tag "${meta_predictions.id}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.75':
        'quay.io/biocontainers/biopython:1.75' }"

    input:
    tuple val(meta_predictions), path(predictions)
    tuple val(meta_spacers), path(spacers_metadata)
    path(metadata)

    output:
    tuple val(meta_predictions), path("custom_chosen_host.tsv"),    emit: tsv
    path "versions.yml",                                      emit: versions

    script:
    """
    collect_crispr_host_info.py \\
        -p ${predictions} \\
        -c ${spacers_metadata} \\
        -m ${metadata} \\
        -o custom_chosen_host.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
