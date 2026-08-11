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
    tuple val(meta), path("${category}.fna"), path("${category}.gff"), path("${category}.faa"),      emit: chosen_sequences
    tuple val("${task.process}"), val('python'), eval('python --version 2>&1 | sed "s/Python //g"'), topic: versions

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
    """
}
