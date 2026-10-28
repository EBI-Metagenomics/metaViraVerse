nextflow.enable.types = true

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

record PlasmidReps {
    meta: Map
    fna: Path
    faa: Path
    gff: Path   // representative GFF: source of the protein-id -> (renamed) contig-id crosswalk
}

record PlasquidResult {
    meta: Map
    fasta: Path                    // combined plasmid contig sequences
    report: Path                   // per-contig RIP/MOB/Inc-group report
    protein_report: Path           // per-protein RIP/MOB/Inc-group report
    rip_seqs_faa: Path             // RIP protein sequences
    mobility_classification: Path  // per-contig conjugative/mobilizable/non_mobilizable
    mobility_stats: Path           // protein_report record counts (JSON): total, RIP_domain/MOB_group/Inc_group present
}

workflow PLASQUID_WORKFLOW {

    take:
    plasmids: Channel<PlasmidReps>

    main:

    //
    // Protein ids keep their original (pre-rename) form, e.g.
    // "MGYG000517684_26|plasmid-1:6000_1", while the matching nucleotide contig has
    // since been renamed to a short accession, e.g. "seq15" -- every downstream
    // script needs this crosswalk to resolve a protein hit back to its contig.
    //
    protein_to_contig = MAP_PROTEIN_TO_CONTIG(
        plasmids.map { r -> record(meta: r.meta, gff: r.gff) }
    )
    reps = plasmids.join(
        protein_to_contig.map { r -> record(meta: r.meta, protein_to_contig: r.map) },
        by: 'meta'
    )

    rep = REPSEARCH(
        reps.map { r -> record(meta: r.meta, proteins: r.faa, protein_to_contig: r.protein_to_contig) },
        file(params.repsearch_db),
        file(params.repfilter_db)
    )

    rna = RNASEARCH(
        plasmids.map { r -> record(meta: r.meta, fna: r.fna) },
        file(params.rna_inc_db)
    )

    inc = INCSEARCH(
        reps.join(rna, by: 'meta')
            .map { r -> record(meta: r.meta, proteins: r.faa, rna_candidates: r.rna_candidates, protein_to_contig: r.protein_to_contig) },
        file(params.incsearch_db)
    )

    mob = MOBSEARCH(
        reps.map { r -> record(meta: r.meta, proteins: r.faa, protein_to_contig: r.protein_to_contig) },
        file(params.mobsearch_db)
    )

    plasquid = EXTRACT_PLASMIDS_DATA(
        inc.join(mob, by: 'meta')
            .join(rep, by: 'meta')
            .join(plasmids, by: 'meta')
            .map { r -> record(meta: r.meta, filt_classification: r.filt_classification, mob_table: r.mob_table, rep_domains: r.rep_domains, fna: r.fna, faa: r.faa) }
    )

    emit:
    plasquid: Channel<PlasquidResult> = plasquid
}
