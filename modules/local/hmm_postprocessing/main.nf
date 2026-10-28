nextflow.enable.types = true

process HMM_POSTPROCESSING {
    /*
    input: File_hmmer_ViPhOG.tbl
    output: File_hmmer_ViPhOG_modified.tbl
    */
    tag "${meta.id}"
    label 'process_single'

    container 'quay.io/microbiome-informatics/virify-python3:1.2'

    input:
    record(
        meta: Map,
        hmmer_tbl: Path
    )

    output:
    record(
        meta: meta,
        modified_tsv: file("${meta.id}_modified.tsv")
    )

    script:
    """
    hmmer_format_table.py -i ${hmmer_tbl} -o ${meta.id}_modified
    """

    stub:
    """
    echo stub > ${meta.id}_modified.tsv
    """
}
