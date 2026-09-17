process EXTRACT_PLASMIDS_DATA {
    /*
     * Combine REPSEARCH/INCSEARCH/MOBSEARCH's independent plasmid-evidence tables
     * into one per-contig report, and extract the reported contigs' own sequences and
     * the protein sequences behind each RIP call.
    */

    label 'process_medium'
    tag "${meta.id}"
    container "docker://mgimenez720/plasquid:latest"

    input:
    tuple val(meta), path(filt_classification), path(mob_table), path(rep_domains), path(fna), path(faa)

    output:
    tuple val(meta), path("plasmids_contigs.fasta"), emit: fasta
    tuple val(meta), path("plasmid_report.tsv"),     emit: report
    tuple val(meta), path("rip_seqs.faa"),           emit: rip_seqs_faa
    path "versions.yml",                             emit: versions

    script:
    """
    plasquid_retrieve_rip_plasmids.R ${filt_classification} ${rep_domains} ${mob_table} ${fna}
    plasquid_rip_extraction.R ${faa} ${filt_classification} ${rep_domains}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | sed -n '1s/^R version \\([0-9.]*\\).*/\\1/p')
    END_VERSIONS
    """
}
