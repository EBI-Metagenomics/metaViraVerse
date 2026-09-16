process MOBSEARCH {
    /*
     * Identify mobilisation/relaxase (MOB) genes: hmmsearch predicted proteins
     * against the MOB profile database and keep hits passing each MOB family's
     * curated score cutoff.
    */

    label 'process_medium'
    tag "${meta.id}"
    container "docker://mgimenez720/plasquid:latest"

    input:
    tuple val(meta), path(proteins)
    path(mobsearch_db)

    output:
    tuple val(meta), path("mob_candidates.tsv"), emit: mob_candidates
    tuple val(meta), path("mob_table.tsv"),      emit: mob_table
    tuple val(meta), path("mob_seqs.faa"),       emit: mob_seqs
    path "versions.yml",                         emit: versions

    script:
    """
    hmmsearch -o log --cpu ${task.cpus} --domtblout mob_candidates.tsv ${mobsearch_db} ${proteins}

    plasquid_filter_mob.R mob_candidates.tsv
    plasquid_mob_extraction.R ${proteins} mob_table.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hmmer: \$(hmmsearch -h | sed -n 's/^# HMMER \\([0-9.]*\\).*/\\1/p')
        r-base: \$(R --version | sed -n '1s/^R version \\([0-9.]*\\).*/\\1/p')
    END_VERSIONS
    """
}
