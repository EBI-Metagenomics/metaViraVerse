process VITAP {

    label 'process_medium'
    tag "$meta.id"

    container "quay.io/microbiome-informatics/vitap:1.7"

    input:
    tuple val(meta), path(fasta)
    path(db)

    output:
    tuple val(meta), path("${meta.id}_vitap_best.tsv"), emit: best_lineages

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

    """
}