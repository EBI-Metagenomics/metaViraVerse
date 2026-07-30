/*
 * modules/local/filter_faa/main.nf
 *
 * Filters .faa protein sequences for ESMFold compatibility:
 *   - Removes sequences longer than max_len residues (default 1000aa)
 *   - Removes sequences with >5% ambiguous residues (X, B, Z, U, O)
 *   - Gracefully handles missing .faa (emits warning, not failure)
 */

process FILTER_FAA {
    tag "${meta.id}"
    label 'process_low'

    conda 'bioconda::biopython=1.83'
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/biopython:1.83' :
        'biocontainers/biopython:1.83' }"

    input:
    tuple val(meta), path(faa)

    output:
    tuple val(meta), path("${meta.id}_filtered.faa"),  emit: filtered_faa
    tuple val(meta), path("${meta.id}_filter_log.tsv"), emit: filter_log

    when:
    task.ext.when == null || task.ext.when

    script:
    def max_len = task.ext.max_len ?: 1000
    """
    if [ ! -s "${faa}" ]; then
        echo "WARNING: No .faa provided or file is empty for ${meta.id} — skipping protein structure subworkflow" >&2
        touch "${meta.id}_filtered.faa"
        printf "seq_id\\tlength\\treason_excluded\\n" > "${meta.id}_filter_log.tsv"
        exit 0
    fi

    python3 ${projectDir}/scripts/filter_faa_for_esmfold.py \\
        --input "${faa}" \\
        --output "${meta.id}_filtered.faa" \\
        --max-len ${max_len} \\
        --log "${meta.id}_filter_log.tsv"
    """

    stub:
    """
    touch "${meta.id}_filtered.faa"
    touch "${meta.id}_filter_log.tsv"
    """
}
