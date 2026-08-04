process MOBSEARCH {

    label 'process_medium'
    tag "${meta.id}"
    container "docker://mgimenez720/plasquid:latest"

    input:
    tuple val(meta), path(proteins)
    path(db_search)

    output:
    tuple val(meta), path("mob_candidates.tsv"), emit: mob_candidates
    tuple val(meta), path("mob_table.tsv"),      emit: mob_table
    tuple val(meta), path("mob_seqs.faa"),       emit: mob_seqs

    script:
    """
    echo "hmsearch"
    hmmsearch -o log --cpu ${task.cpus} --domtblout mob_candidates.tsv ${db_search} ${proteins}

    echo "MOB filter"
    plasquid_filter_mob.R mob_candidates.tsv   # output: mob_table.tsv

    echo "MOB extraction"
    plasquid_mob_extraction.R ${proteins} mob_table.tsv  # output: mob_seqs.faa
    """
}
