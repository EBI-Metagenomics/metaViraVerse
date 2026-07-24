process RNASEARCH {
    label 'process_medium'
    tag "${meta.id}"
    container "docker://mgimenez720/plasquid:latest"

    input:
    tuple val(meta), path(fna)
    path(db_search)

    output:
    tuple val(meta), path("RNA_candidates.tsv"), emit: rna_candidates

    script:
    """

    cmsearch --cpu ${task.cpus} --tblout RNA_candidates.tsv ${db_search} ${fna}

    """
}