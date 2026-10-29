nextflow.enable.types = true

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { KRONA_KTIMPORTTEXT               } from '../../../modules/nf-core/krona/ktimporttext'
include { SANKEY_PLOT                      } from '../../../modules/local/sankey_plot'
include { PLOT_ITOL                        } from '../../../modules/local/plot_itol'


record TaxonomyTables {
    meta: Map
    taxonomy_counts: Path
    metadata_table: Path
}

record TaxonomyPlots {
    meta: Map                 // the input meta, with `tool`
    krona: Path
    sankey: Path?
    itol_tree: Path
    itol_labels: Path?
    itol_counts: Path
    itol_types: Path?
    itol_biomes: Path?
}

workflow TAXONOMY_VISUALISATION {

    take:
    tables: Channel<TaxonomyTables>

    main:

    //
    // -------- Taxonomy visualisation with krona
    //
    krona = KRONA_KTIMPORTTEXT (
        tables.map { r -> tuple(r.meta, r.taxonomy_counts) }
    )

    //
    // -------- Taxonomy visualisation with sankey
    //
    sankey = SANKEY_PLOT (
        tables.map { r -> record(meta: r.meta, table: r.taxonomy_counts) }
    )

    //
    // -------- Taxonomy tree files to upload to iTOL
    //
    itol = PLOT_ITOL (
        tables.map { r -> record(meta: r.meta, table: r.taxonomy_counts, metadata: r.metadata_table) }
    )

    plots = krona.html.map { meta, html -> record(meta: meta, krona: html) }
        .join(sankey, by: 'meta')
        .join(itol, by: 'meta')
        .map { r ->
            record(
                meta: r.meta,
                krona: r.krona,
                sankey: r.sankey_html,
                itol_tree: r.tree_nwk,
                itol_labels: r.tree_labels,
                itol_counts: r.tree_counts,
                itol_types: r.tree_types,
                itol_biomes: r.tree_biomes
            )
        }

    emit:
    plots: Channel<TaxonomyPlots> = plots
}
