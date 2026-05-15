/*
 * subworkflows/local/protein_structure/main.nf
 *
 * Protein structure prediction and annotation subworkflow for metaViraVerse.
 *
 * Triggered only when `.faa` is present in the samplesheet.
 * Existing pipeline behaviour is completely unchanged for users who do not
 * supply protein sequences (caveat 8 — graceful optional handling).
 *
 * Execution order:
 *   FILTER_FAA         → removes seqs >1000aa or >5% ambiguous residues
 *   ESMFOLD            → predicts structures; pLDDT per residue in B-factor
 *   FOLDSEEK_SEARCH    → structural homology vs BFVD (+ optional PDB100)
 *   ECOD_ANNOTATE      → HHsearch-based domain classification from Zenodo DB
 *   MERGE_STRUCT_ANNOT → joins all annotations into extended stats TSV
 *
 * Input:
 *   ch_faa       — channel of [ meta, faa_file ] tuples; may be empty
 *   ch_reps_stats — channel of [ meta, tsv ] from upstream clustering step
 *   bfvd_db      — path or 'BFVD' (streams via bfvd.foldseek.com)
 *   plddt_cutoff — float; minimum mean pLDDT for ECOD/ProteinCartography
 *
 * Output:
 *   struct_stats — extended stats TSV (10 new structural columns)
 *   pdb_dir      — directory of per-protein .pdb files
 *   bfvd_hits    — Foldseek BFVD results TSV
 *   ecod_tsv     — ECOD domain annotation TSV
 */

include { FILTER_FAA          } from '../../../modules/local/filter_faa/main'
include { ESMFOLD             } from '../../../modules/local/esmfold/main'
include { FOLDSEEK_SEARCH     } from '../../../modules/local/foldseek/main'
include { ECOD_ANNOTATE       } from '../../../modules/local/ecod/main'
include { MERGE_STRUCT_ANNOT  } from '../../../modules/local/merge_struct_annotations/main'
include { PROTEINCARTOGRAPHY  } from '../../../modules/local/proteincartography/main'

workflow PROTEIN_STRUCTURE {
    take:
    ch_faa        // channel: [ val(meta), path(faa) ] — may be empty
    ch_reps_stats // channel: [ val(meta), path(tsv) ]
    bfvd_db       // val: path string or 'BFVD'
    plddt_cutoff  // val: float

    main:
    // ── Guard: warn and drop samples with no .faa ────────────
    ch_faa_valid = ch_faa
        .filter { meta, faa ->
            def present = faa && faa.exists() && faa.size() > 0
            if (!present) {
                log.warn "[PROTEIN_STRUCTURE] No .faa for sample '${meta.id}' — " +
                         "skipping protein structure annotation for this sample."
            }
            return present
        }

    // ── Step 1: filter sequences ─────────────────────────────
    FILTER_FAA(ch_faa_valid)

    // ── Step 2: predict structures ───────────────────────────
    ESMFOLD(FILTER_FAA.out.filtered_faa)

    // ── Step 3+4: structural homology search ─────────────────
    FOLDSEEK_SEARCH(
        ESMFOLD.out.pdb_dir,
        bfvd_db
    )

    // ── Step 5: ECOD domain annotation ───────────────────────
    ECOD_ANNOTATE(
        ESMFOLD.out.pdb_dir,
        ESMFOLD.out.confidence_tsv,
        plddt_cutoff
    )

    // ── Step 6: merge all annotations ────────────────────────
    // Join reps_stats with per-sample annotation outputs by meta.id
    ch_merge_input = ch_reps_stats
        .join(ESMFOLD.out.confidence_tsv,      by: [0])
        .join(FOLDSEEK_SEARCH.out.bfvd_hits,   by: [0])
        .join(ECOD_ANNOTATE.out.ecod_tsv,      by: [0])
        .map { meta, reps, conf, bfvd, ecod ->
            [ meta, reps, conf, bfvd, ecod ]
        }

    MERGE_STRUCT_ANNOT(
        ch_merge_input.map { meta, reps, conf, bfvd, ecod -> [ meta, reps  ] },
        ch_merge_input.map { meta, reps, conf, bfvd, ecod -> [ meta, conf  ] },
        ch_merge_input.map { meta, reps, conf, bfvd, ecod -> [ meta, bfvd  ] },
        ch_merge_input.map { meta, reps, conf, bfvd, ecod -> [ meta, ecod  ] }
    )

    // ── Optional: ProteinCartography structural landscape ────
    if (params.run_proteincartography) {
        // taxonomy_tsv: pass empty file if not available
        ch_taxonomy = ch_reps_stats
            .map { meta, tsv -> [ meta, tsv ] }

        PROTEINCARTOGRAPHY(
            ESMFOLD.out.pdb_dir,
            ESMFOLD.out.confidence_tsv,
            ch_taxonomy.map { meta, tsv -> [ meta, tsv ] }
        )
    }

    emit:
    struct_stats = MERGE_STRUCT_ANNOT.out.struct_stats_tsv
    pdb_dir      = ESMFOLD.out.pdb_dir
    bfvd_hits    = FOLDSEEK_SEARCH.out.bfvd_hits
    ecod_tsv     = ECOD_ANNOTATE.out.ecod_tsv
}
