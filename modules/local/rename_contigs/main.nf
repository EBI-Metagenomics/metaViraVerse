/*
 * change headers to short names for contigs fasta and gff
*/
process RENAME_CONTIGS {

    label 'process_low'
    tag "combined"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.75':
        'quay.io/biocontainers/biopython:1.75' }"

    input:
    path(fna)
    path(gff)
    val(start_accession)
    val(end_accession)
    val(types)
    val(biomes)

    output:
    tuple val(["id": "combined"]), path("*renamed*.fasta"),    emit: fna_renamed
    tuple val(["id": "combined"]), path("*renamed*.gff"),      emit: gff_renamed
    tuple val(["id": "combined"]), path("map.tsv"),            emit: map_file
    path "versions.yml",                                     emit: versions

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${params.rename_accession}"
    def fna_args   = fna.collect { it }.join(' ')
    def gff_args   = gff.collect { it }.join(' ')
    def start = start_accession ? "--start ${start_accession}" : ""
    def end = end_accession ? "--end ${end_accession}" : ""
    def type_args  = types.join(' ')
    def biome_args = biomes.join(' ')

    """
    rename_contigs.py \\
       --fasta ${fna_args} \\
       --gff ${gff_args} \\
       --map map.tsv \\
       --prefix ${prefix} \\
       ${start} \\
       ${end} \\
       --type ${type_args} \\
       --biome ${biome_args} \\
       ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        biopython: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('biopython').version)")
    END_VERSIONS
    """
}
