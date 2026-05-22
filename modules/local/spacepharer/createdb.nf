process SPACEPHARER_CREATEDB {

    label 'process_medium'
    tag "$meta.id"

    container "quay.io/biocontainers/spacepharer:5.c2e680a--hd6d6fdc_7"

    input:
    tuple val(meta), path(fasta)
    val(target_db_name)

    output:
    tuple val(meta), path("${target_db_name}*"), emit: db
    path "versions.yml",                        emit: versions

    script:
    def args = task.ext.args   ?: ''
    """
    spacepharer createsetdb \\
      --threads $task.cpus \\
      ${fasta} \\
      ${target_db_name} \\
      tmpFolder \\
      ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spacepharer: \$(spacepharer | grep Version | sed 's/SpacePHARER Version: //g')
    END_VERSIONS
    """
}
