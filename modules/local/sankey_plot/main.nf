process SANKEY_PLOT {

    label 'process_low'
    tag "${meta.id}"
    container 'community.wave.seqera.io/library/pip_plotly:43c7ed2f4992fc92'

    input:
      tuple val(meta), path(table)

    output:
      tuple val(meta), path("*sankey.html"), optional: true, emit: sankey_html
      tuple val("${task.process}"), val('python'), eval("python --version 2>&1 | sed 's/Python //g'"), emit: versions_python, topic: versions
      tuple val("${task.process}"), val('plotly'), eval('python -c "import plotly; print(plotly.__version__)"'), emit: versions_plotly, topic: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    plot_taxonomy_sankey.py \\
        ${args} \\
        --input ${table} \\
        --output ${prefix}_sankey.html
    """
}
