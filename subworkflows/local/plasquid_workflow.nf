/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { REPSEARCH                              } from '../../modules/plasquid/repsearch'
include { RNASEARCH                              } from '../../modules/plasquid/rnasearch'
include { INCSEARCH                              } from '../../modules/plasquid/incsearch'
include { MOBSEARCH                              } from '../../modules/plasquid/mobsearch'
include { EXTRACT_PLASMIDS_DATA                  } from '../../modules/plasquid/extract'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PLASQUID_WORKFLOW {

    take:
    plasmids_fna
    plasmids_faa

    main:

    REPSEARCH(
        plasmids_faa,
        params.repsearch_db,
        params.repfilter_db
    )

    RNASEARCH(
        plasmids_fna,
        params.rna_inc_db
    )

    INCSEARCH(
        plasmids_faa.join(RNASEARCH.out.rna_candidates),
        params.incsearch_db
    )

    MOBSEARCH(
        plasmids_faa,
        params.mobsearch_db
    )

    EXTRACT_PLASMIDS_DATA(
        INCSEARCH.out.filt_classification
           .join(MOBSEARCH.out.mob_table)
           .join(REPSEARCH.out.rep_domains)
           .join(plasmids_fna)
           .join(plasmids_faa)
    )

}

