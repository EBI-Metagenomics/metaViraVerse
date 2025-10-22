process SANKEY_PLOT {

    label 'process_low'
    tag "${meta.id} ${set_name}"
    container 'quay.io/microbiome-informatics/sankeyd3:0.12.3'

    input:
      tuple val(meta), path(table)

    output:
      tuple val(meta), val(set_name), path("*.sankey.html")

    script:

    """
    sankey_plot.R ${table} "${meta.id}.${set_name}.sankey.html"
    """
}