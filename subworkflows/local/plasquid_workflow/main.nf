/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MAP_PROTEIN_TO_CONTIG                  } from '../../../modules/plasquid/map_protein_to_contig'
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
    plasmids_gff   // representative GFF: source of the protein-id -> (renamed) contig-id crosswalk

    main:

    ch_versions = channel.empty()

    //
    // Protein ids keep their original (pre-rename) form, e.g.
    // "MGYG000517684_26|plasmid-1:6000_1", while the matching nucleotide contig has
    // since been renamed to a short accession, e.g. "seq15" -- every downstream
    // script needs this crosswalk to resolve a protein hit back to its contig.
    //
    MAP_PROTEIN_TO_CONTIG(
        plasmids_gff
    )
    ch_versions = ch_versions.mix(MAP_PROTEIN_TO_CONTIG.out.versions)

    REPSEARCH(
        plasmids_faa.join(MAP_PROTEIN_TO_CONTIG.out.map),
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
        plasmids_faa
           .join(RNASEARCH.out.rna_candidates)
           .join(MAP_PROTEIN_TO_CONTIG.out.map),
        params.incsearch_db
    )
    ch_versions = ch_versions.mix(INCSEARCH.out.versions)

    MOBSEARCH(
        plasmids_faa.join(MAP_PROTEIN_TO_CONTIG.out.map),
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
    fasta          = EXTRACT_PLASMIDS_DATA.out.fasta           // combined plasmid contig sequences
    report         = EXTRACT_PLASMIDS_DATA.out.report           // per-contig RIP/MOB/Inc-group report
    protein_report = EXTRACT_PLASMIDS_DATA.out.protein_report   // per-protein RIP/MOB/Inc-group report
    rip_seqs       = EXTRACT_PLASMIDS_DATA.out.rip_seqs_faa      // RIP protein sequences
    mobility_stats = EXTRACT_PLASMIDS_DATA.out.mobility_stats   // conjugative/mobilizable/non_mobilizable counts (JSON)
    versions       = ch_versions
}
