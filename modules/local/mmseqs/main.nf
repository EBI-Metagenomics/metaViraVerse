process MMSEQS {

    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fe/fe49c17754753d6cd9a31e5894117edaf1c81e3d6053a12bf6dc8f3af1dffe23/data'
        : 'community.wave.seqera.io/library/mmseqs2:18.8cc5c--af05c9a98d9f6139'}"

    input:
    tuple val(meta), path(faa_file)
    val id_threshold
    val cov_threshold

    output:
    path "protein_catalogue-*.tsv", emit: mmseq_cluster_tsv
    tuple val("${task.process}"), val('mmseqs'), eval('mmseqs version'), topic: versions, emit: versions_mmseqs

    script:
    int threshold_rounded = id_threshold * 100;
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
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch protein_catalogue-90.tsv
    """
}