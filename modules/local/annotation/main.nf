nextflow.enable.types = true

process ANNOTATION {
    /*
     * Generate tabular file with ViPhOG annotation results for
     * proteins predicted in viral contigs
    */
    tag "${meta.id}"

    label 'process_single'

    container 'quay.io/microbiome-informatics/virify-python3:1.1'

    input:
    record(
        meta: Map,
        tab: Path,
        gff: Path
    )

    output:
    record(
        meta: meta,
        annotations: file("*_annotation.tsv")
    )

    script:
    """
    viral_contigs_annotation.py -o . -t ${tab} -g ${gff}
    """

    stub:
    """
    echo stub > stub_annotation.tsv
    """
}
