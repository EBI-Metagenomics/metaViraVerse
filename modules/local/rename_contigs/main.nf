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
    val(meta)
    path(fna)
    path(gff)
    val(types)
    val(biomes)
    val(start_accession)
    val(end_accession)
    path(fna_tp)
    path(gff_tp)
    val(types_tp)
    val(biomes_tp)
    val(sources_tp)
    val(viral_sequence_identifier)
    val(prophage_identifier)
    val(plasmid_identifier)

    output:
    tuple val(meta), path("combined.fna"), emit: fna_renamed
    tuple val(meta), path("combined.gff"), emit: gff_renamed
    tuple val(meta), path("combined.tsv"), emit: map_file
    tuple val("${task.process}"), val('python'), eval('python --version 2>&1 | sed "s/Python //g"'), topic: versions

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${params.rename_accession}"

    def fna_files    = fna    ? fna.collect { it }    : []
    def gff_files    = gff    ? gff.collect { it }    : []
    def fna_tp_files = fna_tp ? fna_tp.collect { it } : []
    def gff_tp_files = gff_tp ? gff_tp.collect { it } : []

    def fna_args    = fna_files    ? "--fasta ${fna_files.join(' ')}"       : ''
    def gff_args    = gff_files    ? "--gff ${gff_files.join(' ')}"         : ''
    def fna_tp_args = fna_tp_files ? "--fasta-tp ${fna_tp_files.join(' ')}" : ''
    def gff_tp_args = gff_tp_files ? "--gff-tp ${gff_tp_files.join(' ')}"   : ''

    def start = start_accession ? "--start ${start_accession}" : ''
    def end   = end_accession   ? "--end ${end_accession}"     : ''

    def type_args      = types      ? "--type ${types.join(' ')}"           : ''
    def biome_args     = biomes     ? "--biome ${biomes.join(' ')}"         : ''
    def type_tp_args   = types_tp   ? "--type-tp ${types_tp.join(' ')}"     : ''
    def biome_tp_args  = biomes_tp  ? "--biome-tp ${biomes_tp.join(' ')}"   : ''
    def source_tp_args = sources_tp ? "--source-tp ${sources_tp.join(' ')}" : ''

    def viral_id_arg    = viral_sequence_identifier ? "--viral-sequence-identifier ${viral_sequence_identifier}" : ''
    def prophage_id_arg = prophage_identifier       ? "--prophage-identifier ${prophage_identifier}"             : ''
    def plasmid_id_arg  = plasmid_identifier        ? "--plasmid-identifier ${plasmid_identifier}"               : ''

    """
    rename_contigs.py \\
       ${fna_args} \\
       ${gff_args} \\
       ${fna_tp_args} \\
       ${gff_tp_args} \\
       --map combined.tsv \\
       --prefix ${prefix} \\
       ${start} \\
       ${end} \\
       ${type_args} \\
       ${biome_args} \\
       ${type_tp_args} \\
       ${biome_tp_args} \\
       ${source_tp_args} \\
       ${viral_id_arg} \\
       ${prophage_id_arg} \\
       ${plasmid_id_arg} \\
       ${args}
    """
}
