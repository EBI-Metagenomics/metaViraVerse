nextflow.enable.types = true

process BUILD_FINAL_GFF {

    label 'process_low'
    tag "${meta.id}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.75':
        'quay.io/biocontainers/biopython:1.75' }"

    input:
    record(
        meta: Map,
        input_gff: Path
    )
    bacphlip: Path
    hmmfiles: List<Path>    // hmmer tables, empty when --annotation_db is not set
    additional_gff: Path?   // AMR GFF, null when all AMR tools are skipped

    output:
    record(
        meta: meta,
        final_gff: file("*final.gff*")
    )

    topic:
    tuple(task.process, 'python', eval("python --version 2>&1 | sed 's/Python //g'")) >> 'versions'

    script:
    // hmmer tables and AMR GFF are optional upstream (empty when the step was skipped)
    def hmmer_arg = hmmfiles       ? "--hmmer ${hmmfiles.join(' ')}" : ""
    def gff_arg   = additional_gff ? "--gff ${additional_gff}"       : ""
    """
    build_final_gff.py \\
        --input ${input_gff} \\
        --output ${meta.id}_final.gff \\
        --bacphlip ${bacphlip} \\
        ${hmmer_arg} \\
        ${gff_arg} \\
        --compress-output
    """
}
