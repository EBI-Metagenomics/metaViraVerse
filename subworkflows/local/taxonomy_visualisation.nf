/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { IQTREE                           } from '../../modules/nf-core/iqtree'
include { KRONA_KTIMPORTTEXT               } from '../../modules/nf-core/krona/ktimporttext'
include { MAFFT_ALIGN                      } from '../../modules/nf-core/mafft/align'

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

    if ( params.iqtree ) {
        //
        // -------- Alignment for IQTree
        //
        MAFFT_ALIGN(
            fna_reps_seqs,
            fna_reps_seqs.map { meta, fasta -> tuple([:], []) },
            fna_reps_seqs.map { meta, fasta -> tuple([:], []) },
            fna_reps_seqs.map { meta, fasta -> tuple([:], []) },
            fna_reps_seqs.map { meta, fasta -> tuple([:], []) },
            fna_reps_seqs.map { meta, fasta -> tuple([:], []) },
            true
        )

        //
        // -------- IQTree
        //
        IQTREE(
            MAFFT_ALIGN.out.fas.map { meta, align -> tuple(meta, align, []) },
            [], [], [], [], [], [], [], [], [], [], [], []
        )
        ch_versions = ch_versions.mix(IQTREE.out.versions)
    }

    emit:

    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}