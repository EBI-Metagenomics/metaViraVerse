nextflow.enable.types = true

process SPACEPHARER_EASYPREDICT {

    label 'process_medium'
    tag "$meta.id"

    container "quay.io/biocontainers/spacepharer:5.c2e680a--hd6d6fdc_7"

    input:
    record(
        meta: Map,
        crispr_fasta: Path,
        target_db: Set<Path>,
        target_db_rev: Set<Path>
    )
    target_db_name: String

    output:
    record(
        meta: meta,
        predictions: file("${meta.id}_predictions.tsv")
    )

    topic:
    tuple(task.process, 'spacepharer', eval("spacepharer | grep Version | sed 's/SpacePHARER Version: //g'")) >> 'versions'

    script:
    """
    spacepharer easy-predict \\
      --threads $task.cpus \\
      ${crispr_fasta} \\
      ${target_db_name} \\
      ${meta.id}_predictions.tsv \\
      tmp
    """

    stub:
    """
    echo stub > ${meta.id}_predictions.tsv
    """
}
