process EXTRACT {

    label 'process_medium'
    tag "${meta.id}"
    container "docker://mgimenez720/plasquid:latest"

    input:
    tuple val(meta), path(filt_tsv), path(domains), path(mob_table), path(fna), path(faa)

    output:
    tuple val(meta), path("Plasmids_contigs.fasta"), emit: fasta
    tuple val(meta), path("Plasmid_Report.tsv"), emit: report
    tuple val(meta), path("RIP_seqs.faa"), emit: rip_seqs_faa
    tuple val(meta), path("Result.tsv"), emit: repsearch_result
    tuple val(meta), path("Result.fasta"), emit: repsearch_fasta

    script:
    """
    echo "GeneRetrieve"
    Retrieve_RIP_plasmids.R ${filt_tsv} ${domains} ${mob_table} ${fna}

    echo "RipExtract"
    RIP_extraction.R ${faa} ${filt_tsv} ${domains}

    echo "RepSearchOut"
    repsearch_sum.R Plasmid_Report.tsv ${fna}
    """
}