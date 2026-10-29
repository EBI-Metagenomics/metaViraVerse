nextflow.enable.types = true

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
    record(
        meta: Map,
        reps_file: Path
    )
    full_gff: Path
    mapfile: Path?

    output:
    record(
        meta: meta,
        reps_list: file("${meta.id}_reps.tsv"),
        reps_gff: file("${meta.id}_reps.gff"),
        reps_proteins_list: file("${meta.id}_reps_proteins.tsv")
    )

    topic:
    tuple(task.process, 'python', eval("python --version 2>&1 | sed 's/Python //g'")) >> 'versions'
    tuple(task.process, 'biopython', eval("python -c \"import pkg_resources; print(pkg_resources.get_distribution('biopython').version)\"")) >> 'versions'

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

    stub:
    """
    echo stub > ${meta.id}_reps.tsv
    echo stub > ${meta.id}_reps.gff
    echo stub > ${meta.id}_reps_proteins.tsv
    """
}
