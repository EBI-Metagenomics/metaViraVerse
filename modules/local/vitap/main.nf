nextflow.enable.types = true

process VITAP {

    label 'process_medium'
    tag "$meta.id"

    container "quay.io/microbiome-informatics/vitap:1.7"

    input:
    record(
        meta: Map,
        fasta: Path
    )
    db: Path

    output:
    record(
        meta: meta,
        best_lineages: file("*_vitap_best.tsv")
    )

    // version is hardcoded because Docker container has version 1.7 installed
    // but tool itself still prints version 1.3 (which is a bug in tool)
    topic:
    tuple(task.process, 'VITAP', '1.7') >> 'versions'

    script:
    def fasta_file = fasta.name.endsWith('.gz') ? fasta.baseName : fasta.name
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    if [[ ${fasta} == *.gz ]]; then
        gunzip -c ${fasta} > ${fasta_file}
    fi

    VITAP assignment \\
      -i ${fasta_file} \\
      -d ${db} \\
      -o ${prefix}_vitap

    # filter records (because VITAP adds random genomic fragments for normalization and calibration)
    # take input sequence names
    grep '>' ${fasta_file} | sed 's/>//' > names.txt

    # write header first
    head -n1 ${prefix}_vitap/best_determined_lineages.tsv > ${prefix}_vitap_best.tsv

    # add filtered records
    grep -w -f names.txt ${prefix}_vitap/best_determined_lineages.tsv >> ${prefix}_vitap_best.tsv

    # Attention:
    # version is hardcoded because Docker container has version 1.7 installed
    # but tool itself still prints version 1.3 (which is a bug in tool)
    """
}
