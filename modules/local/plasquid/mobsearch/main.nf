nextflow.enable.types = true

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
    record(
        meta: Map,
        proteins: Path,
        protein_to_contig: Path
    )
    mobsearch_db: Path

    output:
    record(
        meta: meta,
        mob_candidates: file("mob_candidates.tsv"),
        mob_table: file("mob_table.tsv"),
        mob_seqs: file("mob_seqs.faa")
    )

    topic:
    tuple(task.process, 'hmmer', eval("hmmsearch -h | sed -n 's/^# HMMER \\([0-9.]*\\).*/\\1/p'")) >> 'versions'
    tuple(task.process, 'r-base', eval("R --version | sed -n '1s/^R version \\([0-9.]*\\).*/\\1/p'")) >> 'versions'

    script:
    """
    hmmsearch -o log --cpu ${task.cpus} --domtblout mob_candidates.tsv ${mobsearch_db} ${proteins}

    plasquid_filter_mob.R mob_candidates.tsv ${protein_to_contig}
    plasquid_mob_extraction.R ${proteins} mob_table.tsv
    """
}
