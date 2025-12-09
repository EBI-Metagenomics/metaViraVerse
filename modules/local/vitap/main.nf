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
    """
    gunzip -c "${fasta}" > input.fasta

    VITAP assignment \\
      -i input.fasta \\
      -d ${db} \\
      -o ${meta.id}_vitap

    # filter records
    grep '>' input.fasta | sed 's/>//' > names.txt

    # header
    head -n1 ${meta.id}_vitap/best_determined_lineages.tsv > ${meta.id}_vitap_best.tsv

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