/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MAP_PROTEIN_TO_CONTIG                  } from '../../../modules/local/plasquid/map_protein_to_contig'
include { REPSEARCH                              } from '../../../modules/local/plasquid/repsearch'
include { RNASEARCH                              } from '../../../modules/local/plasquid/rnasearch'
include { INCSEARCH                              } from '../../../modules/local/plasquid/incsearch'
include { MOBSEARCH                              } from '../../../modules/local/plasquid/mobsearch'
include { EXTRACT_PLASMIDS_DATA                  } from '../../../modules/local/plasquid/extract_plasmids_data'
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

    //
    // Protein ids keep their original (pre-rename) form, e.g.
    // "MGYG000517684_26|plasmid-1:6000_1", while the matching nucleotide contig has
    // since been renamed to a short accession, e.g. "seq15" -- every downstream
    // script needs this crosswalk to resolve a protein hit back to its contig.
    //
    MAP_PROTEIN_TO_CONTIG(
        plasmids_gff.map { meta, gff -> record(meta: meta, gff: gff) }
    )
    protein_to_contig = MAP_PROTEIN_TO_CONTIG.out.map { r -> [r.meta, r.map] }

    REPSEARCH(
        plasmids_faa
           .join(protein_to_contig)
           .map { meta, proteins, map -> record(meta: meta, proteins: proteins, protein_to_contig: map) },
        params.repsearch_db,
        params.repfilter_db
    )

    RNASEARCH(
        plasmids_fna.map { meta, fna -> record(meta: meta, fna: fna) },
        params.rna_inc_db
    )

    INCSEARCH(
        plasmids_faa
           .join(RNASEARCH.out.map { r -> [r.meta, r.rna_candidates] })
           .join(protein_to_contig)
           .map { meta, proteins, rna, map -> record(meta: meta, proteins: proteins, rna_candidates: rna, protein_to_contig: map) },
        params.incsearch_db
    )

    MOBSEARCH(
        plasmids_faa
           .join(protein_to_contig)
           .map { meta, proteins, map -> record(meta: meta, proteins: proteins, protein_to_contig: map) },
        params.mobsearch_db
    )

    EXTRACT_PLASMIDS_DATA(
        INCSEARCH.out.map { r -> [r.meta, r.filt_classification] }
           .join(MOBSEARCH.out.map { r -> [r.meta, r.mob_table] })
           .join(REPSEARCH.out.map { r -> [r.meta, r.rep_domains] })
           .join(plasmids_fna)
           .join(plasmids_faa)
           .map { meta, filt, mob, rep, fna, faa -> record(meta: meta, filt_classification: filt, mob_table: mob, rep_domains: rep, fna: fna, faa: faa) }
    )
    plasmids_data = EXTRACT_PLASMIDS_DATA.out

    emit:
    fasta                  = plasmids_data.map { r -> [r.meta, r.fasta] }                   // combined plasmid contig sequences
    report                 = plasmids_data.map { r -> [r.meta, r.report] }                   // per-contig RIP/MOB/Inc-group report
    protein_report         = plasmids_data.map { r -> [r.meta, r.protein_report] }           // per-protein RIP/MOB/Inc-group report
    rip_seqs               = plasmids_data.map { r -> [r.meta, r.rip_seqs_faa] }              // RIP protein sequences
    mobility_classification = plasmids_data.map { r -> [r.meta, r.mobility_classification] } // per-contig conjugative/mobilizable/non_mobilizable
    mobility_stats         = plasmids_data.map { r -> [r.meta, r.mobility_stats] }           // protein_report record counts (JSON): total, RIP_domain/MOB_group/Inc_group present
}
