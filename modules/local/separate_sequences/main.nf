process SEPARATE_SEQUENCES {

    label 'process_low'
    tag "$meta.id"

    input:
    tuple val(meta), path(fasta)
    val pattern

    output:
    tuple val(meta), path("${meta.id}_${pattern}.fa"), emit: chosen_sequences

    script:
    """
    awk '/^>/{f=(\$0 ~ /${pattern}/)} f' ${fasta} > ${meta.id}_${pattern}.fa
    """
}
