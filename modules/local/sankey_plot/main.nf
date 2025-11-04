process SANKEY_PLOT {

    label 'process_low'
    tag "${meta.id}"
    container 'community.wave.seqera.io/library/pip_plotly:43c7ed2f4992fc92'

    input:
      tuple val(meta), path(table)

    output:
      tuple val(meta), path("*sankey.html"), optional: true, emit: sankey_html
      path "versions.yml",                                   emit: versions

    script:

    """
    plot_taxonomy_sankey.py \\
        --input ${table} \\
        --output ${meta.id}_sankey.html \\
        --input-format tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        plotly: \$(python -c "import plotly; print(plotly.__version__)")
    END_VERSIONS
    """
}