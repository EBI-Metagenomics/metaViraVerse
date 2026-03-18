/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { KRONA_KTIMPORTTEXT               } from '../../modules/nf-core/krona/ktimporttext'
include { SANKEY_PLOT                      } from '../../modules/local/sankey_plot'


workflow TAXONOMY_VISUALISATION {

    take:
    fna_reps_seqs
    reps_krona_tsv

    main:

    ch_versions = channel.empty()

    //
    // -------- Taxonomy visualisation with krona
    //
    KRONA_KTIMPORTTEXT (
        reps_krona_tsv
    )
    ch_versions = ch_versions.mix(KRONA_KTIMPORTTEXT.out.versions)

    //
    // -------- Taxonomy visualisation with sankey
    //
    SANKEY_PLOT (
        reps_krona_tsv
    )
    ch_versions = ch_versions.mix(SANKEY_PLOT.out.versions)

    emit:

    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}