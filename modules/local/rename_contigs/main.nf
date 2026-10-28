nextflow.enable.types = true

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
    meta: Map
    fna: List<Path>
    gff: List<Path>
    types: List<String>
    biomes: List<String>
    start_accession: Integer?
    end_accession: Integer?
    fna_tp: List<Path>
    gff_tp: List<Path>
    types_tp: List<String>
    biomes_tp: List<String>
    sources_tp: List<String>
    viral_sequence_identifier: String?
    prophage_identifier: String?
    plasmid_identifier: String?
    rename_accession: String

    output:
    record(
        meta: meta,
        fna_renamed: file("combined.fna"),
        gff_renamed: file("combined.gff"),
        map_file: file("combined.tsv")
    )

    topic:
    tuple(task.process, 'python', eval('python --version 2>&1 | sed "s/Python //g"')) >> 'versions'

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${rename_accession}"

    def fna_files    = fna
    def gff_files    = gff
    def fna_tp_files = fna_tp
    def gff_tp_files = gff_tp

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

    stub:
    """
    echo stub > combined.fna
    echo stub > combined.gff
    echo stub > combined.tsv
    """
}
