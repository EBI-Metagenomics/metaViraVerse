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
    tuple val(meta), path(proteins), path(protein_to_contig)
    path(repsearch_db)
    path(repfilter_db)

    output:
    tuple val(meta), path("prots_vs_rep.tsv"), emit: prots_vs_rep
    tuple val(meta), path("rep_domains.tsv"),  emit: rep_domains
    path "versions.yml",                       emit: versions

    script:
    """
    hmmsearch --cut_ga --cpu ${task.cpus} -o log --domtblout prots_vs_rep.tsv ${repsearch_db} ${proteins}

    plasquid_dom_arch.R prots_vs_rep.tsv
    plasquid_filter_rip.R single_dom_rip.tsv domain_architecture.RDS ${repfilter_db} ${protein_to_contig}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hmmer: \$(hmmsearch -h | sed -n 's/^# HMMER \\([0-9.]*\\).*/\\1/p')
        r-base: \$(R --version | sed -n '1s/^R version \\([0-9.]*\\).*/\\1/p')
    END_VERSIONS
    """
}
