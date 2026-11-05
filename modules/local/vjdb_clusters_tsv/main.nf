nextflow.enable.types = true

/*
 * Convert the vjdb-clustering membership CSV (vjdb_id, copy_of, ..., cluster_rep) to the
 * vclust cluster TSV (object, cluster) that the downstream steps read.
 *
 * CHOOSE_SEQUENCES has already removed identical sequences (sha256 of the upper-cased
 * sequence), which is stricter than the case-sensitive md5 vjdb-clustering uses. So vjdb
 * must find no copies: every row has copy_of == vjdb_id. A copy here means the two hashes
 * disagree, and a later update (which reads seq_hash from this CSV) would see it too, so
 * the step fails instead of passing it on.
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

    output:
    record(
        meta: meta,
        clusters_tsv: file("${meta.id}_clusters.tsv"),
        clusters_csv: clusters_csv,
        reps_fasta: reps_fasta
    )

    script:
    """
    awk -F',' -v OFS='\\t' '
        NR == 1 {
            for (i = 1; i <= NF; i++) col[\$i] = i
            if (!("vjdb_id" in col) || !("copy_of" in col) || !("cluster_rep" in col)) {
                print "ERROR: ${clusters_csv} needs the columns vjdb_id, copy_of and cluster_rep" > "/dev/stderr"
                exit 1
            }
            print "object", "cluster"
            next
        }
        \$col["copy_of"] != \$col["vjdb_id"] {
            copies++
            if (copies <= 10) examples = examples " " \$col["vjdb_id"] "=" \$col["copy_of"]
        }
        { print \$col["vjdb_id"], \$col["cluster_rep"] }
        END {
            if (copies > 0) {
                print "ERROR: vjdb-clustering found " copies " identical sequence(s) that CHOOSE_SEQUENCES kept, e.g." examples > "/dev/stderr"
                exit 1
            }
        }
    ' ${clusters_csv} > ${meta.id}_clusters.tsv
    """

    stub:
    """
    printf "object\\tcluster\\n" > ${meta.id}_clusters.tsv
    """
}
