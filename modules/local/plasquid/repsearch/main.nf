nextflow.enable.types = true

process REPSEARCH {
    /*
     * Identify plasmid replication initiator proteins (RIPs): hmmsearch predicted
     * proteins against the RIP profile database, resolve each candidate's
     * single-/multi-domain architecture, and apply curated per-domain score
     * (and domain-architecture) cutoffs to call the real RIPs.
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
    repsearch_db: Path
    repfilter_db: Path

    output:
    record(
        meta: meta,
        prots_vs_rep: file("prots_vs_rep.tsv"),
        rep_domains: file("rep_domains.tsv")
    )

    topic:
    tuple(task.process, 'hmmer', eval("hmmsearch -h | sed -n 's/^# HMMER \\([0-9.]*\\).*/\\1/p'")) >> 'versions'
    tuple(task.process, 'r-base', eval("R --version | sed -n '1s/^R version \\([0-9.]*\\).*/\\1/p'")) >> 'versions'

    script:
    """
    hmmsearch --cut_ga --cpu ${task.cpus} -o log --domtblout prots_vs_rep.tsv ${repsearch_db} ${proteins}

    plasquid_dom_arch.R prots_vs_rep.tsv
    plasquid_filter_rip.R single_dom_rip.tsv domain_architecture.RDS ${repfilter_db} ${protein_to_contig}
    """

    stub:
    """
    echo stub > prots_vs_rep.tsv
    echo stub > rep_domains.tsv
    """
}
