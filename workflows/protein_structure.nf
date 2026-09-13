/*
 * workflows/protein_structure.nf
 *
 * Top-level entry point for the protein structure subworkflow.
 * Called from main.nf when params.run_protein_structure == true
 * and .faa is present in the samplesheet.
 *
 * Add to main.nf:
 *
 *   include { PROTEIN_STRUCTURE } from './workflows/protein_structure'
 *
 *   // After existing samplesheet parsing, extract optional .faa channel:
 *   ch_faa = ch_samplesheet
 *       .map { meta, reads, faa -> faa ? [ meta, faa ] : null }
 *       .filter { it != null }
 *
 *   // Invoke subworkflow if requested:
 *   if (params.run_protein_structure) {
 *       PROTEIN_STRUCTURE(
 *           ch_faa,
 *           params.bfvd_db      ?: 'BFVD',
 *           params.plddt_cutoff ?: 0.7,
 *           VIRALCLUSTER.out.reps_stats
 *       )
 *   }
 *
 * New params to add to nextflow.config / nextflow_schema.json:
 *   run_protein_structure  (boolean, default: false)
 *   esmfold_mode           (string: 'api' | 'local', default: 'api')
 *   bfvd_db                (string: path or 'BFVD', default: 'BFVD')
 *   plddt_cutoff           (number, default: 0.7)
 *   foldseek_search_pdb    (boolean, default: false)
 */

include { PROTEIN_STRUCTURE as PROTEIN_STRUCTURE_WF } from '../subworkflows/local/protein_structure/main'

workflow PROTEIN_STRUCTURE {
    take:
    ch_faa
    bfvd_db
    plddt_cutoff
    reps_stats

    main:
    PROTEIN_STRUCTURE_WF(
        ch_faa,
        bfvd_db,
        plddt_cutoff,
        reps_stats
    )

    emit:
    struct_stats = PROTEIN_STRUCTURE_WF.struct_stats
    pdb_dir      = PROTEIN_STRUCTURE_WF.pdb_dir
    bfvd_hits    = PROTEIN_STRUCTURE_WF.bfvd_hits
}
