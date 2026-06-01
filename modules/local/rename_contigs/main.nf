/*
 * change headers to short names for contigs fasta and gff
*/
process RENAME_CONTIGS {

    label 'process_low'
    tag "${meta.id}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.75':
        'quay.io/biocontainers/biopython:1.75' }"

    input:
    tuple val(meta), path(fna), path(gff)
    val(start_accession)
    val(end_accession)

    output:
    tuple val(meta), path("*_renamed.fasta"),    emit: contigs_renamed
    tuple val(meta), path("*_renamed.gff"),      emit: gff_renamed, optional: true
    tuple val(meta), path("${meta.id}.map.tsv"), emit: map_file
    path "versions.yml",                         emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def gff_arg = gff ? "--gff ${gff}" ? ""
    def start = start_accession ? "--start {start_accession}" ? ""
    def end = end_accession ? "--end {end_accession}" ? ""
    """
    rename_contigs.py \\
       --input ${fna} \\
       ${gff_arg} \\
       --map ${meta.id}.map.tsv \\
       --prefix ${prefix} \\
       ${start} \\
       ${end}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        biopython: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('biopython').version)")
    END_VERSIONS
    """
}
