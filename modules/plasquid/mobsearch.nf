process MOBSEARCH {

    label 'process_medium'
    tag "${meta.id}"
    container "docker://mgimenez720/plasquid:latest"

    input:
    tuple val(meta), path(proteins)
    path(db_search)

    output:
    tuple val(meta), path("Mob_candidates.tsv"), emit: mob_candidates
    tuple val(meta), path("Mob_table.tsv"), emit: mob_table
    tuple val(meta), path("MOB_seqs.faa"), emit: mob_seqs

    script:
    """
    echo "hmsearch"
    hmmsearch -o log --cpu ${task.cpus} --domtblout Mob_candidates.tsv ${db_search} ${proteins}

    echo "MOB filder"
    Filter_Mob.R Mob_candidates.tsv

    echo "MOB extraction"
    MOB_extraction.R ${proteins} Mob_table.tsv
    """
}
