nextflow.enable.types = true

/*
 * Convert the vjdb-clustering membership CSV (vjdb_id, ..., cluster_rep) to the
 * vclust cluster TSV (object, cluster) that the downstream steps read
*/
process VJDB_CLUSTERS_TSV {

    label 'process_single'
    tag "${meta.id}"
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.75':
        'quay.io/biocontainers/biopython:1.75' }"

    input:
    meta: Map
    clusters_csv: Path

    output:
    record(
        meta: meta,
        clusters_tsv: file("${meta.id}_clusters.tsv"),
        clusters_csv: clusters_csv
    )

    script:
    """
    awk -F',' -v OFS='\\t' '
        NR == 1 { for (i = 1; i <= NF; i++) col[\$i] = i; print "object", "cluster"; next }
        { print \$col["vjdb_id"], \$col["cluster_rep"] }
    ' ${clusters_csv} > ${meta.id}_clusters.tsv
    """

    stub:
    """
    printf "object\\tcluster\\n" > ${meta.id}_clusters.tsv
    """
}
