nextflow.enable.types = true

process PLOT_PROTEIN_CLUSTERS_STATS {
    /*
     * Plots and statistics for protein clustering from PhaMMseqs or MMseqs2
     * (bin/plot_protein_clusters_stats.py): cluster size distribution, rank-abundance, top functions,
     * function purity per cluster, and clustering agreement metrics (MI, NMI, AMI, ARI,
     * homogeneity, completeness, V-measure) between clusters and protein functions
     * (or an optional reference labelling).
    */

    label 'process_low'
    tag "${meta.id}"

    conda "${moduleDir}/environment.yml"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/matplotlib_numpy_pandas_python_scikit-learn:9699a1c5d7878366':
        'community.wave.seqera.io/library/matplotlib_numpy_pandas_python_scikit-learn:50a467f1b0184adf' }"

    input:
    record(
        meta: Map,
        clusters: Path,
        fasta: Path?,       // optional: protein FASTA with functions in headers (>proteinID function)
        reference: Path?    // optional: proteinID<TAB>label TSV to compare clusters against; null to use protein functions
    )
    cluster_format: String  // 'phammseqs' (phamID<TAB>proteinID function) or 'mmseqs' (representative<TAB>member)

    output:
    record(
        meta: meta,
        size_distribution_tsv: file("*_cluster_size_distribution.tsv"),
        size_distribution_png: file("*_cluster_size_distribution.png"),
        rank_abundance: file("*_cluster_rank_abundance.png"),
        top_functions_tsv: file("*_top_functions.tsv"),
        top_functions_png: file("*_top_functions.png", optional: true),
        purity_tsv: file("*_cluster_function_purity.tsv"),
        purity_png: file("*_cluster_function_purity.png", optional: true),
        stats: file("*_clustering_stats.tsv")
    )

    topic:
    tuple(task.process, 'python', eval("python --version 2>&1 | sed 's/Python //g'")) >> 'versions'
    tuple(task.process, 'matplotlib', eval('python -c "import matplotlib; print(matplotlib.__version__)"')) >> 'versions'
    tuple(task.process, 'scikit-learn', eval('python -c "import sklearn; print(sklearn.__version__)"')) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference_arg = reference ? "--reference ${reference}" : ""
    def fasta_arg     = fasta     ? "--fasta ${fasta}"         : ""
    if (!(cluster_format in ['phammseqs', 'mmseqs'])) {
        error "PLOT_PROTEIN_CLUSTERS_STATS: cluster_format must be 'phammseqs' or 'mmseqs', got '${cluster_format}'"
    }
    """
    plot_protein_clusters_stats.py \\
        --${cluster_format} ${clusters} \\
        ${fasta_arg} \\
        -o . \\
        --prefix ${prefix} \\
        ${reference_arg} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_cluster_size_distribution.tsv ${prefix}_cluster_size_distribution.png
    touch ${prefix}_cluster_rank_abundance.png
    touch ${prefix}_top_functions.tsv ${prefix}_top_functions.png
    touch ${prefix}_cluster_function_purity.tsv ${prefix}_cluster_function_purity.png
    touch ${prefix}_clustering_stats.tsv
    """
}
