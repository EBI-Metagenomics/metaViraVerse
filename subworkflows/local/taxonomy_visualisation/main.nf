/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { KRONA_KTIMPORTTEXT               } from '../../../modules/nf-core/krona/ktimporttext'
include { SANKEY_PLOT                      } from '../../../modules/local/sankey_plot'
include { PLOT_ITOL                        } from '../../../modules/local/plot_itol'


workflow TAXONOMY_VISUALISATION {

    take:
    reps_krona_and_metadata

    main:

    //
    // -------- Taxonomy visualisation with krona
    //
    KRONA_KTIMPORTTEXT (
        reps_krona_and_metadata.map{meta, tsv, _metadata -> [meta, tsv]}
    )

    //
    // -------- Taxonomy visualisation with sankey
    //
    SANKEY_PLOT (
        reps_krona_and_metadata.map{meta, tsv, _metadata -> [meta, tsv]}
    )

    //
    // -------- Taxonomy tree files to upload to iTOL
    //
    PLOT_ITOL (
        reps_krona_and_metadata
    )
}
