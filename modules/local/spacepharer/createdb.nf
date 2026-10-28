nextflow.enable.types = true

process SPACEPHARER_CREATEDB {

    label 'process_medium'
    tag "$meta.id"

    container "quay.io/biocontainers/spacepharer:5.c2e680a--hd6d6fdc_7"

    input:
    record(
        meta: Map,
        fasta: Path
    )
    target_db_name: String

    output:
    record(
        meta: meta,
        db: files("${target_db_name}*")
    )

    topic:
    tuple(task.process, 'spacepharer', eval("spacepharer | grep Version | sed 's/SpacePHARER Version: //g'")) >> 'versions'

    script:
    def args = task.ext.args   ?: ''
    """
    spacepharer createsetdb \\
      --threads $task.cpus \\
      ${fasta} \\
      ${target_db_name} \\
      tmpFolder \\
      ${args}
    """
}
