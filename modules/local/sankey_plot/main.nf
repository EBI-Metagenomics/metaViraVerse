nextflow.enable.types = true

process SANKEY_PLOT {

    label 'process_low'
    tag "${meta.id}"
    container 'community.wave.seqera.io/library/pip_plotly:43c7ed2f4992fc92'

    input:
    record(
        meta: Map,
        table: Path
    )

    output:
    record(
        meta: meta,
        sankey_html: file("*sankey.html", optional: true),
        sankey_tsv: file("*sankey.tsv", optional: true)     // the links the plot is drawn from
    )

    topic:
    tuple(task.process, 'python', eval("python --version 2>&1 | sed 's/Python //g'")) >> 'versions'
    tuple(task.process, 'plotly', eval('python -c "import plotly; print(plotly.__version__)"')) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    plot_taxonomy_sankey.py \\
        ${args} \\
        --input ${table} \\
        --output ${prefix}_sankey.html \\
        --output-table ${prefix}_sankey.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo stub > ${prefix}_sankey.html
    printf "source_level\\tsource\\ttarget_level\\ttarget\\ttarget_lineage\\tsequences\\n" > ${prefix}_sankey.tsv
    """
}
