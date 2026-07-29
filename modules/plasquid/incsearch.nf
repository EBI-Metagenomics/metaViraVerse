process INCSEARCH {

    label 'process_medium'
    tag "${meta.id}"
    container "docker://mgimenez720/plasquid:latest"

    input:
    tuple val(meta), path(proteins), path(rna_candidates)
    path(db_search)

    output:
    tuple val(meta), path("Inc_candidates.tsv"), emit: inc_candidates
    tuple val(meta), path("Classification_table.tsv"), emit: inc_classification
    tuple val(meta), path("Filtered_Classif.tsv"), emit: filt_classification

    script:
    """
    echo "hmsearch"
    hmmsearch -o log --cpu ${task.cpus} --domtblout Inc_candidates.tsv ${db_search} ${proteins}

    echo "Inc_classification"
    Inc_classification.R Inc_candidates.tsv ${rna_candidates}

    echo "Filter_classification"
    Filter_classification.R Classification_table.tsv
    """
}
