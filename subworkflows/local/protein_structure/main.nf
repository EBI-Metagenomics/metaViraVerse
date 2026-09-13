/*
 * subworkflows/local/protein_structure/main.nf
 *
 * Protein structure annotation subworkflow.
 * Accepts either:
 *   - Pre-computed PDB files (e.g. from nf-core/proteinfold) via ch_pdbs
 *   - Protein FASTA via ch_faa, in which case ESMFold runs first
 *
 * Downstream: Foldseek vs BFVD, ECOD domain annotation,
 * ProteinCartography landscape, merged stats.
 */

include { FILTER_FAA            } from '../../../modules/local/filter_faa/main'
include { ESMFOLD               } from '../../../modules/local/esmfold/main'
include { FOLDSEEK_SEARCH       } from '../../../modules/local/foldseek/main'
include { ECOD_ANNOTATE         } from '../../../modules/local/ecod/main'
include { MERGE_STRUCT_ANNOT    } from '../../../modules/local/merge_struct_annotations/main'
include { PROTEINCARTOGRAPHY    } from '../../../modules/local/proteincartography/main'

workflow PROTEIN_STRUCTURE {

    take:
    ch_faa       // tuple val(meta), path(faa)  — protein FASTA (optional if ch_pdbs provided)
    ch_pdbs      // tuple val(meta), path(pdbs) — pre-computed PDB dir (optional if ch_faa provided)

    main:

    // ── Stage 1: Get PDB structures ──────────────────────────
    // If pre-computed PDBs provided, use them directly
    // Otherwise, predict with ESMFold from FAA

    if (ch_pdbs) {
        ch_structures = ch_pdbs
    } else {
        FILTER_FAA ( ch_faa )
        ESMFOLD ( FILTER_FAA.out.filtered_faa )
        ch_structures = ESMFOLD.out.pdbs
    }

    ch_confidence = ch_pdbs
        ? Channel.empty()
        : ESMFOLD.out.confidence_tsv

    // ── Stage 2: Structural homology search ──────────────────
    FOLDSEEK_SEARCH ( ch_structures )

    // ── Stage 3: Evolutionary domain annotation ──────────────
    ECOD_ANNOTATE ( ch_structures )

    // ── Stage 4: Merge annotations ───────────────────────────
    ch_merge_input = ch_structures
        .join(ch_confidence.ifEmpty { [[],[]] },  by: [0], remainder: true)
        .join(FOLDSEEK_SEARCH.out.bfvd_tsv,       by: [0])
        .join(ECOD_ANNOTATE.out.ecod_tsv,          by: [0])
        .map { meta, structs, conf, bfvd, ecod ->
            [ meta, structs, conf ?: [], bfvd, ecod ]
        }

    MERGE_STRUCT_ANNOT (
        ch_merge_input.map { meta, structs, conf, bfvd, ecod -> [ meta, structs ] },
        ch_merge_input.map { meta, structs, conf, bfvd, ecod -> [ meta, conf    ] },
        ch_merge_input.map { meta, structs, conf, bfvd, ecod -> [ meta, bfvd    ] },
        ch_merge_input.map { meta, structs, conf, bfvd, ecod -> [ meta, ecod    ] }
    )

    // ── Stage 5: Structural landscape (optional) ─────────────
    if (params.run_proteincartography) {
        PROTEINCARTOGRAPHY ( ch_structures )
    }

    emit:
    struct_stats  = MERGE_STRUCT_ANNOT.out.merged_tsv
    bfvd_hits     = FOLDSEEK_SEARCH.out.bfvd_tsv
    ecod_tsv      = ECOD_ANNOTATE.out.ecod_tsv
    structures    = ch_structures
}
