nextflow.enable.types = true

process SORT_GFF {

    label 'process_low'
    tag "$meta.id"

    input:
    record(
        meta: Map,
        gff: Path
    )

    output:
    record(
        meta: meta,
        sorted_gff: file("*sorted*")
    )

    script:
    """
    sort -k1,1V -k4,4n ${gff} > ${gff.baseName}.sorted.gff
    """
}
