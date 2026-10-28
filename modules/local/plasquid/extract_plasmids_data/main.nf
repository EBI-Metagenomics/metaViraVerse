nextflow.enable.types = true

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
    record(
        meta: Map,
        filt_classification: Path,
        mob_table: Path,
        rep_domains: Path,
        fna: Path,
        faa: Path
    )

    output:
    record(
        meta: meta,
        fasta: file("plasmids_contigs.fasta"),
        report: file("plasmid_report.tsv"),
        protein_report: file("protein_report.tsv"),
        rip_seqs_faa: file("rip_seqs.faa"),
        mobility_classification: file("mobility_classification.tsv"),
        mobility_stats: file("mobility_stats.json")
    )

    topic:
    tuple(task.process, 'r-base', eval("R --version | sed -n '1s/^R version \\([0-9.]*\\).*/\\1/p'")) >> 'versions'

    script:
    """
    # Output: plasmids_contigs.fasta, plasmid_report.tsv
    plasquid_retrieve_rip_plasmids.R ${filt_classification} ${rep_domains} ${mob_table} ${fna}

    # Output: rip_seqs.faa, protein_report.tsv, mobility_classification.tsv, mobility_stats.json
    # (mobility_stats.json counts protein_report.tsv records: total, and with RIP_domain/MOB_group/Inc_group)
    plasquid_rip_extraction.R ${faa} ${filt_classification} ${rep_domains} ${mob_table}
    """
}
