process SORT_GFF {

    label 'process_low'
    tag "$meta.id"

    input:
    tuple val(meta), path(gff)

    output:
    tuple val(meta), path("*sorted*"), emit: sorted_gff

    script:
    """
    sort -k1,1V -k4,4n ${gff} > ${gff.baseName}.sorted.gff
    """
}
