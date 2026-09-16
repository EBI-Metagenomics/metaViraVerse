/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { REPSEARCH                              } from '../../../modules/plasquid/repsearch'
include { RNASEARCH                              } from '../../../modules/plasquid/rnasearch'
include { INCSEARCH                              } from '../../../modules/plasquid/incsearch'
include { MOBSEARCH                              } from '../../../modules/plasquid/mobsearch'
include { EXTRACT_PLASMIDS_DATA                  } from '../../../modules/plasquid/extract_plasmids_data'
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

    ch_versions = channel.empty()

    REPSEARCH(
        plasmids_faa,
        params.repsearch_db,
        params.repfilter_db
    )
    ch_versions = ch_versions.mix(REPSEARCH.out.versions)

    RNASEARCH(
        plasmids_fna,
        params.rna_inc_db
    )
    ch_versions = ch_versions.mix(RNASEARCH.out.versions)

    INCSEARCH(
        plasmids_faa.join(RNASEARCH.out.rna_candidates),
        params.incsearch_db
    )
    ch_versions = ch_versions.mix(INCSEARCH.out.versions)

    MOBSEARCH(
        plasmids_faa,
        params.mobsearch_db
    )
    ch_versions = ch_versions.mix(MOBSEARCH.out.versions)

    EXTRACT_PLASMIDS_DATA(
        INCSEARCH.out.filt_classification
           .join(MOBSEARCH.out.mob_table)
           .join(REPSEARCH.out.rep_domains)
           .join(plasmids_fna)
           .join(plasmids_faa)
    )
    ch_versions = ch_versions.mix(EXTRACT_PLASMIDS_DATA.out.versions)

    emit:
    fasta      = EXTRACT_PLASMIDS_DATA.out.fasta        // combined plasmid contig sequences
    report     = EXTRACT_PLASMIDS_DATA.out.report        // per-contig RIP/MOB/Inc-group report
    rip_seqs   = EXTRACT_PLASMIDS_DATA.out.rip_seqs_faa   // RIP protein sequences
    versions   = ch_versions
}
