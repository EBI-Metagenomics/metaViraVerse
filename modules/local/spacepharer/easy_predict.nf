process SPACEPHARER_EASYPREDICT {

    label 'process_medium'
    tag "$meta.id"

    container "quay.io/biocontainers/spacepharer:5.c2e680a--hd6d6fdc_7"

    input:
    tuple val(meta), path(crispr_fasta)
    path(target_db)
    path(target_db_rev)

    output:
    tuple val(meta), path("${meta.id}_predictions.tsv"), emit: predictions
    path "versions.yml",                                 emit: versions

    script:
    """
    spacepharer easy-predict \\
      --threads $task.cpus \\
      ${crispr_fasta} \\
      ${target_db} \\
      ${meta.id}_predictions.tsv \\
      tmp

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spacepharer: \$(spacepharer | grep Version | sed 's/SpacePHARER Version: //g')
    END_VERSIONS
    """
}
