process VITAP {

    label 'process_medium'
    tag "$meta.id"

    container "quay.io/microbiome-informatics/vitap:1.7"

    input:
    tuple val(meta), path(fasta)
    path(db)

    output:
    tuple val(meta), path("${meta.id}_vitap_best.tsv"), emit: best_lineages
    path "versions.yml",                                emit: versions

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
    head -n1 ${meta.id}_vitap/best_determined_lineages.tsv > ${meta.id}_vitap_best.tsv

    # add filtered records
    grep -w -f names.txt ${meta.id}_vitap/best_determined_lineages.tsv >> ${meta.id}_vitap_best.tsv

    # Attention:
    # version is hardcoded because Docker container has version 1.7 installed
    # but tool itself still prints version 1.3 (which is a bug in tool)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        VITAP: 1.7
    END_VERSIONS
    """
}
