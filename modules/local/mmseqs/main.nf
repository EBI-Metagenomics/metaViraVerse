nextflow.enable.types = true

process MMSEQS {

    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fe/fe49c17754753d6cd9a31e5894117edaf1c81e3d6053a12bf6dc8f3af1dffe23/data'
        : 'community.wave.seqera.io/library/mmseqs2:18.8cc5c--af05c9a98d9f6139'}"

    input:
    record(
        meta: Map,
        faa_file: Path
    )
    id_threshold: Float
    cov_threshold: Float

    output:
    record(
        meta: meta,
        mmseq_cluster_tsv: file("protein_catalogue-*.tsv")
    )

    topic:
    tuple(task.process, 'mmseqs', eval('mmseqs version')) >> 'versions'

    script:
    def threshold_rounded = Math.round(id_threshold * 100)
    """
    timestamp() {
        date +"%H:%M:%S"
    }
    echo "\$(timestamp) [mmseqs script] Creating MMseqs database"

    mmseqs createdb ${faa_file} mmseqs.db

    echo "\$(timestamp) [mmseqs script] Clustering MMseqs with linclust with option -c ${id_threshold}"

    mmseqs linclust \
    mmseqs.db \
    mmseqs_cluster.db \
    mmseqs-tmp --min-seq-id ${id_threshold} \
    --threads ${task.cpus} \
    -c ${cov_threshold} \
    --cov-mode 1 \
    --cluster-mode 2 \
    --kmer-per-seq 80

    echo "\$(timestamp) [mmseqs script] Parsing output to create FASTA file of all sequences"

    mmseqs createseqfiledb mmseqs.db \
    mmseqs_cluster.db \
    mmseqs_cluster_seq \
    --threads ${task.cpus}

    mmseqs result2flat mmseqs.db \
    mmseqs.db \
    mmseqs_cluster_seq \
    mmseqs_cluster.fa

    echo "\$(timestamp) [mmseqs script] Parsing output to create TSV file with cluster membership"

    mmseqs createtsv mmseqs.db \
    mmseqs.db \
    mmseqs_cluster.db \
    protein_catalogue-${threshold_rounded}.tsv \
    --threads ${task.cpus}
    """

    stub:
    """
    touch protein_catalogue-90.tsv
    """
}
