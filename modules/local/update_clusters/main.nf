nextflow.enable.types = true

/*
 * --mode update with an existing catalogue or an update that is not clustered with
 * vjdb-clustering (bin/update_clusters.py).
*/

/*
 * The existing catalogue was clustered with --viral_clustering default, the update uses
 * vjdb: build the vjdb membership CSV and representative FASTA of the existing clusters.
*/
process CLUSTERS_TO_VJDB {

    label 'process_single'
    tag "${meta.id}"
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.75':
        'quay.io/biocontainers/biopython:1.75' }"

    input:
    record(
        meta: Map,
        clusters_tsv: Path,   // clusters of the existing catalogue
        fasta: Set<Path>      // every viral sequence of the existing catalogue (viral sequences, prophages)
    )

    stage:
    stageAs fasta, 'sequences/*'

    output:
    record(
        meta: meta,
        clusters_csv: file("${meta.id}_existing_vjdb_clusters.csv"),
        reps_fasta: file("${meta.id}_existing_vjdb_reps.fna.gz")
    )

    script:
    """
    update_clusters.py to-vjdb \\
        --clusters-tsv ${clusters_tsv} \\
        --fasta ${fasta.join(' ')} \\
        --prefix ${meta.id}_existing
    """

    stub:
    """
    echo "vjdb_id,copy_of,seq_gc,seq_len,seq_hash,cluster_rep" > ${meta.id}_existing_vjdb_clusters.csv
    echo "" | gzip > ${meta.id}_existing_vjdb_reps.fna.gz
    """
}

/*
 * The update clusters the whole catalogue again (--viral_clustering default): the clusters
 * whose representative was not one before are annotated, the existing representatives that
 * no longer are are dropped from the merged annotations.
*/
process COMPARE_CLUSTERS {

    label 'process_single'
    tag "${meta.id}"
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.75':
        'quay.io/biocontainers/biopython:1.75' }"

    input:
    record(
        meta: Map,
        new_clusters: Path    // clusters of the whole updated catalogue
    )
    old_clusters: Path        // clusters of the existing catalogue

    stage:
    stageAs old_clusters, 'existing/*'

    output:
    record(
        meta: meta,
        annotate_tsv: file("${meta.id}_clusters_to_annotate.tsv"),
        dropped_reps: file("${meta.id}_dropped_reps.txt"),
        cluster_events: file("${meta.id}_cluster_events.tsv")
    )

    script:
    """
    update_clusters.py compare \\
        --old-clusters ${old_clusters} \\
        --new-clusters ${new_clusters} \\
        --prefix ${meta.id}
    """

    stub:
    """
    printf "object\\tcluster\\n" > ${meta.id}_clusters_to_annotate.tsv
    touch ${meta.id}_dropped_reps.txt
    printf "cluster_rep\\tevent\\told_cluster_reps\\tn_existing\\tn_new\\n" > ${meta.id}_cluster_events.tsv
    """
}
