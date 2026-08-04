process EXTRACT_PLASMIDS_DATA {

    label 'process_medium'
    tag "${meta.id}"
    container "docker://mgimenez720/plasquid:latest"

    input:
    tuple val(meta), path(filt_tsv), path(mob_table), path(domains), path(fna), path(faa)

    output:
    tuple val(meta), path("plasmids_contigs.fasta"), emit: fasta
    tuple val(meta), path("plasmid_report.tsv"),     emit: report
    tuple val(meta), path("rip_seqs.faa"),           emit: rip_seqs_faa
    tuple val(meta), path("result.tsv"),             emit: repsearch_result
    tuple val(meta), path("result.fasta"),           emit: repsearch_fasta

    script:
    """
    echo "GeneRetrieve step"
    plasquid_retrieve_rip_plasmids.R ${filt_tsv} ${domains} ${mob_table} ${fna}  # output:  plasmids_contigs.fasta, plasmid_report.tsv

    echo "RipExtract step"
    plasquid_rip_extraction.R ${faa} ${filt_tsv} ${domains}  # output: rip_seqs.faa

    echo "RepSearchOut step"
    plasquid_repsearch_sum.R plasmid_report.tsv ${fna}  # output: result.tsv, result.fasta
    """
}
