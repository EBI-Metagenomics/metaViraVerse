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
    tuple val(meta), path(proteins), path(protein_to_contig)
    path(mobsearch_db)

    output:
    tuple val(meta), path("mob_candidates.tsv"), emit: mob_candidates
    tuple val(meta), path("mob_table.tsv"),      emit: mob_table
    tuple val(meta), path("mob_seqs.faa"),       emit: mob_seqs
    tuple val("${task.process}"), val('hmmer'), eval("hmmsearch -h | sed -n 's/^# HMMER \\([0-9.]*\\).*/\\1/p'"), emit: versions_hmmer, topic: versions
    tuple val("${task.process}"), val('r-base'), eval("R --version | sed -n '1s/^R version \\([0-9.]*\\).*/\\1/p'"), emit: versions_r_base, topic: versions

    script:
    """
    hmmsearch -o log --cpu ${task.cpus} --domtblout mob_candidates.tsv ${mobsearch_db} ${proteins}

    plasquid_filter_mob.R mob_candidates.tsv ${protein_to_contig}
    plasquid_mob_extraction.R ${proteins} mob_table.tsv
    """
}
