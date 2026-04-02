process BUILD_FINAL_GFF {

    label 'process_low'
    tag "combined"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.75':
        'quay.io/biocontainers/biopython:1.75' }"

    input:
    tuple val(meta_input), path(input_gff)
    tuple val(meta_bacphlip), path(bacphlip)
    tuple val(meta_hmmer), path(hmmfile)
    tuple val(meta_gff), path(additional_gff)

    output:
    path("*.gff"),             emit: final_gff
    path "versions.yml",       emit: versions

    script:
    """
    build_final_gff.py \\
        --input ${input_gff} \\
        --output ${meta_input.id}_final.gff \\
        --bacphlip ${bacphlip} \\
        --hmmer ${hmmfile} \\
        --gff ${additional_gff} \\
        --compress-output

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
