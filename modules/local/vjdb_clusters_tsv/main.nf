nextflow.enable.types = true

/*
 * Convert the vjdb-clustering membership CSV (vjdb_id, copy_of, ..., cluster_rep) to the
 * vclust cluster TSV (object, cluster) that the downstream steps read (bin/vjdb_clusters_tsv.py).
 *
 * With the cluster events of a vjdb update, also lists the clusters to annotate (only the
 * new clusters: grown and merged clusters keep an old, already annotated representative)
 * and the old representatives dropped by merges. Without events every cluster is annotated.
 *
 * Fails when vjdb found identical sequences: CHOOSE_SEQUENCES must have removed them.
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
    reps_fasta: Path
    events: Path?          // cluster_events.tsv of a vjdb update; null for the initial clustering

    output:
    record(
        meta: meta,
        clusters_tsv: file("${meta.id}_clusters.tsv"),
        annotate_tsv: file("${meta.id}_clusters_to_annotate.tsv"),
        dropped_reps: file("${meta.id}_dropped_reps.txt"),
        clusters_csv: clusters_csv,
        reps_fasta: reps_fasta
    )

    script:
    def events_arg = events ? "--events ${events}" : ''
    """
    vjdb_clusters_tsv.py \\
        --clusters-csv ${clusters_csv} \\
        --prefix ${meta.id} \\
        ${events_arg}
    """

    stub:
    """
    printf "object\\tcluster\\n" > ${meta.id}_clusters.tsv
    printf "object\\tcluster\\n" > ${meta.id}_clusters_to_annotate.tsv
    touch ${meta.id}_dropped_reps.txt
    """
}
