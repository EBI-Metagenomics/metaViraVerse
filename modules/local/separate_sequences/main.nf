nextflow.enable.types = true

process SEPARATE_SEQUENCES {

    label 'process_low'
    tag "${meta.id}_${category}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.75':
        'quay.io/biocontainers/biopython:1.75' }"

    input:
    record(
        meta: Map,
        fna: Path
    )
    gff: Path?
    faa: Path?
    map_file: Path
    category: String

    output:
    record(
        meta: meta,
        fna: file("${category}.fna"),
        gff: file("${category}.gff"),
        faa: file("${category}.faa")
    )

    topic:
    tuple(task.process, 'python', eval('python --version 2>&1 | sed "s/Python //g"')) >> 'versions'

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

    stub:
    """
    echo stub > ${category}.fna
    echo stub > ${category}.gff
    echo stub > ${category}.faa
    """
}
