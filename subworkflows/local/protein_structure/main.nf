/*
 * subworkflows/local/protein_structure/main.nf
 *
 * Protein structure prediction and annotation subworkflow for metaViraVerse.
 * Triggered only when `.faa` is present in the samplesheet.
 * Existing pipeline behaviour is unchanged for users without protein sequences.
 *
 * Input:
 *   ch_faa       — channel of [ meta, faa_file ] tuples (optional; may be empty)
 *   bfvd_db      — path to Foldseek BFVD database directory (or 'BFVD' to stream)
 *   plddt_cutoff — float: minimum mean pLDDT for ECOD/ProteinCartography (default 0.7)
 *   reps_stats   — path to existing viral_sequences_reps_stats.tsv
 *
 * Output:
 *   struct_stats — enriched stats TSV with 9 structural columns
 *   pdb_dir      — directory of predicted .pdb files
 *   bfvd_hits    — Foldseek BFVD results TSV
 */

include { FILTER_FAA          } from '../../../modules/local/filter_faa/main'
include { ESMFOLD             } from '../../../modules/local/esmfold/main'
include { FOLDSEEK_SEARCH     } from '../../../modules/local/foldseek/main'
include { ECOD_ANNOTATE       } from '../../../modules/local/ecod/main'
include { MERGE_STRUCT_ANNOT  } from '../../../modules/local/merge_struct_annotations/main'

workflow PROTEIN_STRUCTURE {
    take:
    ch_faa           // channel: [ val(meta), path(faa) ] — may be empty Channel
    bfvd_db          // path
    plddt_cutoff     // val(float)
    reps_stats       // channel: [ val(meta), path(tsv) ]

    main:
    // Guard: warn and short-circuit if no .faa provided
    ch_faa_valid = ch_faa.filter { meta, faa ->
        def has_faa = faa && faa.exists() && faa.size() > 0
        if (!has_faa) {
            log.warn "[PROTEIN_STRUCTURE] No .faa provided for sample ${meta.id} — skipping protein structure subworkflow"
        }
        return has_faa
    }

    // Step 1 — filter sequences for ESMFold compatibility
    FILTER_FAA(ch_faa_valid)

    // Step 2 + 3 — predict structures + implicit pLDDT QC in confidence TSV
    ESMFOLD(FILTER_FAA.out.filtered_faa)

    // Step 4 + 5 — structural homology search vs BFVD (+ optional PDB)
    FOLDSEEK_SEARCH(
        ESMFOLD.out.pdb_dir,
        bfvd_db
    )

    // Step 6 + 7 — ECOD domain annotation + SCOP cross-refs
    ECOD_ANNOTATE(
        ESMFOLD.out.pdb_dir,
        ESMFOLD.out.confidence_tsv,
        plddt_cutoff
    )

    // Step 8 — merge all annotations into extended stats TSV
    MERGE_STRUCT_ANNOT(
        reps_stats,
        ESMFOLD.out.confidence_tsv,
        FOLDSEEK_SEARCH.out.bfvd_hits,
        ECOD_ANNOTATE.out.ecod_tsv
    )

    emit:
    struct_stats = MERGE_STRUCT_ANNOT.out.struct_stats_tsv
    pdb_dir      = ESMFOLD.out.pdb_dir
    bfvd_hits    = FOLDSEEK_SEARCH.out.bfvd_hits
    ecod_tsv     = ECOD_ANNOTATE.out.ecod_tsv
}
