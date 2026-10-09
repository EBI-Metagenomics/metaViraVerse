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
    tuple val(meta), path(clusters)
    val(cluster_format) // 'phammseqs' (phamID<TAB>proteinID function) or 'mmseqs' (representative<TAB>member)
    path(fasta)         // optional: protein FASTA with functions in headers (>proteinID function); [] to skip
    path(reference)     // optional: proteinID<TAB>label TSV to compare clusters against; [] to use protein functions

    output:
    tuple val(meta), path("${prefix}_cluster_size_distribution.tsv"), path("${prefix}_cluster_size_distribution.png"), emit: size_distribution
    tuple val(meta), path("${prefix}_cluster_rank_abundance.png"),                                                     emit: rank_abundance
    tuple val(meta), path("${prefix}_top_functions.tsv"),                                                              emit: top_functions_tsv
    tuple val(meta), path("${prefix}_top_functions.png"), optional: true,                                              emit: top_functions_png
    tuple val(meta), path("${prefix}_cluster_function_purity.tsv"),                                                    emit: purity_tsv
    tuple val(meta), path("${prefix}_cluster_function_purity.png"), optional: true,                                    emit: purity_png
    tuple val(meta), path("${prefix}_clustering_stats.tsv"),                                                       emit: stats
    path "versions.yml",                                                                                           emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        matplotlib: \$(python -c "import matplotlib; print(matplotlib.__version__)")
        scikit-learn: \$(python -c "import sklearn; print(sklearn.__version__)")
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_cluster_size_distribution.tsv ${prefix}_cluster_size_distribution.png
    touch ${prefix}_cluster_rank_abundance.png
    touch ${prefix}_top_functions.tsv ${prefix}_top_functions.png
    touch ${prefix}_cluster_function_purity.tsv ${prefix}_cluster_function_purity.png
    touch ${prefix}_clustering_stats.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        matplotlib: \$(python -c "import matplotlib; print(matplotlib.__version__)")
        scikit-learn: \$(python -c "import sklearn; print(sklearn.__version__)")
    END_VERSIONS
    """
}
