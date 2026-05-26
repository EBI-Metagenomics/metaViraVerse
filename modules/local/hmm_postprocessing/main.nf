process HMM_POSTPROCESSING {
    /*
    input: File_hmmer_ViPhOG.tbl
    output: File_hmmer_ViPhOG_modified.tbl
    */
    tag "${meta.id}"
    label 'process_single'
    
    container 'quay.io/microbiome-informatics/virify-python3:1.2'

    input:
      tuple val(meta), path(hmmer_tbl)
    
    output:
      tuple val(meta), path("${meta.id}_modified.tsv")
    
    script:
    """
    hmmer_format_table.py -i ${hmmer_tbl} -o ${meta.id}_modified -t ${params.hmmer_tool}
    """
}


