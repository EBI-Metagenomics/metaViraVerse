/*
 * module to extruct information about cluster reps from initial GFF
*/
process EXTRACT_REPS_STATS {

    label 'process_low'
    tag "${meta.id}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.75':
        'quay.io/biocontainers/biopython:1.75' }"

    input:
    tuple val(meta), path(reps_file)
    path(full_gff)
    path(mapfile)

    output:
    tuple val(meta), path("${meta.id}_reps.tsv"),           emit: reps_list
    tuple val(meta), path("${meta.id}_reps.gff"),           emit: reps_gff
    tuple val(meta), path("${meta.id}_reps_proteins.tsv"),  emit: reps_proteins_list
    tuple val("${task.process}"), val('python'), eval("python --version 2>&1 | sed 's/Python //g'"), emit: versions_python, topic: versions
    tuple val("${task.process}"), val('biopython'), eval("python -c \"import pkg_resources; print(pkg_resources.get_distribution('biopython').version)\""), emit: versions_biopython, topic: versions

    script:
    def mapfile_arg = mapfile ? "--mapfile ${mapfile}" : ''
    """
    extract_reps_stats.py \\
       --viral-list ${reps_file} \\
       --gff ${full_gff} \\
       --output-reps-list ${meta.id}_reps.tsv \\
       --output-reps-gff ${meta.id}_reps.gff \\
       --output-reps-proteins ${meta.id}_reps_proteins.tsv \\
       ${mapfile_arg}
    """
}
