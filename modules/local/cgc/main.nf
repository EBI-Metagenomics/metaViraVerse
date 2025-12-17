process COMBINED_GENE_CALLER {

    label 'process_low'
    tag "${meta.id}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.75':
        'quay.io/biocontainers/biopython:1.75' }"

    input:
    tuple val(meta), path(prodigal_faa)
    tuple val(meta), path(phanotate_faa)

    output:
    tuple val(meta), path("*_cgc.faa"),    emit: combined_proteins
    path "versions.yml",                   emit: versions

    script:
    """
    combined_gene_caller.py \\
       -p ${prodigal_faa} \\
       -t ${phanotate_faa} \\
       -o ${meta.id}_cgc.faa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        biopython: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('biopython').version)")
    END_VERSIONS
    """
}
