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
    tuple val(meta), path(full_gff), path(mapfile)
    tuple val(meta2), path(reps_file)

    output:
    tuple val(meta), path("${meta.id}_reps_stats.tsv"), emit: reps_stats_tsv
    tuple val(meta), path("${meta.id}_krona.tsv"),      emit: reps_krona_tsv
    path "versions.yml",                                emit: versions

    script:
    """
    extract_reps_stats.py \\
       --viral-list ${reps_file} \\
       --gff ${full_gff} \\
       --mapfile ${mapfile} \\
       --output ${meta.id}_reps_stats.tsv \\
       --krona ${meta.id}_krona.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        biopython: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('biopython').version)")
    END_VERSIONS
    """
}