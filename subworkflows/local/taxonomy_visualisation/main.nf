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

workflow TAXONOMY_VISUALISATION {

    take:
    tables: Channel<TaxonomyTables>

    main:

    //
    // -------- Taxonomy visualisation with krona
    //
    KRONA_KTIMPORTTEXT (
        tables.map { r -> tuple(r.meta, r.taxonomy_counts) }
    )

    //
    // -------- Taxonomy visualisation with sankey
    //
    SANKEY_PLOT (
        tables.map { r -> record(meta: r.meta, table: r.taxonomy_counts) }
    )

    //
    // -------- Taxonomy tree files to upload to iTOL
    //
    PLOT_ITOL (
        tables.map { r -> record(meta: r.meta, table: r.taxonomy_counts, metadata: r.metadata_table) }
    )
}
