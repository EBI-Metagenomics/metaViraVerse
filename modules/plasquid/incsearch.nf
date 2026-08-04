process INCSEARCH {

    label 'process_medium'
    tag "${meta.id}"
    container "docker://mgimenez720/plasquid:latest"

    input:
    tuple val(meta), path(proteins), path(rna_candidates)
    path(db_search)

    output:
    tuple val(meta), path("inc_candidates.tsv"),       emit: inc_candidates
    tuple val(meta), path("classification_table.tsv"), emit: inc_classification
    tuple val(meta), path("filtered_classif.tsv"),     emit: filt_classification

    script:
    """
    echo "hmsearch"
    hmmsearch -o log --cpu ${task.cpus} --domtblout inc_candidates.tsv ${db_search} ${proteins}

    echo "Inc_classification"
    plasquid_inc_classification.R inc_candidates.tsv ${rna_candidates}

    echo "Filter_classification"
    plasquid_filter_classification.R classification_table.tsv
    """
}
