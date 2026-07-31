process SEPARATE_SEQUENCES {

    label 'process_low'
    tag "${meta.id}_${category}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.75':
        'quay.io/biocontainers/biopython:1.75' }"


    input:
    tuple val(meta), path(fna)
    tuple val(meta_gff), path(gff)
    tuple val(meta_faa), path(faa)
    tuple val(meta_map), path(map_file)
    val category

    output:
    tuple val(meta), path("${category}.fna"),                emit: chosen_sequences
    tuple val(meta), path("${category}.gff"), optional: true, emit: chosen_gff
    tuple val(meta), path("${category}.faa"), optional: true, emit: chosen_faa
    path "versions.yml",                                      emit: versions

    script:
    def gff_arg = gff ? "--gff ${gff}" : ""
    def faa_arg = faa ? "--faa ${faa}" : ""
    """
    separate_sequences.py \\
       --fna ${fna} \\
       ${gff_arg} \\
       ${faa_arg} \\
       --map ${map_file} \\
       --category ${category} \\
       --output-prefix ${category}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        biopython: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('biopython').version)")
    END_VERSIONS
    """
}
