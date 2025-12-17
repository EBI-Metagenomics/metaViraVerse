// https://github.com/chg60/phamclust

process PHAMCLUST {

    label 'process_medium'
    tag "${meta.id}"
    container "docker: quay.io/microbiome-informatics/phamclust:latest"

    input:
    tuple val(meta), path(table)

    output:
    tuple val(meta), path("${meta.id}_phamclust"),    emit: phamclust_results
    path "versions.yml",                              emit: versions

    script:
    """
    phamclust ${table} ${meta.id}_phamclust

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phamclust: commit e032f49
    END_VERSIONS
    """
}