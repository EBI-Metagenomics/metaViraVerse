nextflow.enable.types = true

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { GUNZIP as UNCOMPRESSED_REPS_FNA              } from '../../../modules/nf-core/gunzip'
include { GUNZIP as UNCOMPRESSED_REPS_FAA              } from '../../../modules/nf-core/gunzip'
include { SEQTK_SUBSEQ as GREP_FNA                     } from '../../../modules/nf-core/seqtk/subseq'
include { SEQTK_SUBSEQ as GREP_FAA                     } from '../../../modules/nf-core/seqtk/subseq'

include { EXTRACT_REPS_STATS                           } from '../../../modules/local/extract_reps_stats'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

record ClusterRepsInput {
    meta: Map
    clusters_tsv: Path
}

record SequencesInput {
    meta: Map
    fna: Path
}

record ClusterReps {
    meta: Map
    reps_list: Path            // representative ids
    reps_gff: Path
    reps_proteins_list: Path
    fna: Path                  // representative sequences, uncompressed
    faa: Path                  // representative proteins, uncompressed
    fna_gz: Path
    faa_gz: Path
}

workflow EXTRACT_CLUSTER_FILES {

    take:
    clusters: Channel<ClusterRepsInput>
    sequences: Channel<SequencesInput>
    gff: Value<Path>
    faa: Value<Path>
    mapfile: Value<Path>?   // rename map (null for plasmids)
    type: String

    main:

    //
    // ----------- Statistics and taxonomy from GFF for reps -----------
    //

    reps_stats = EXTRACT_REPS_STATS (
        clusters.map { r -> record(meta: r.meta, reps_file: r.clusters_tsv) },
        gff,
        mapfile
    )

    //
    // ----------- Extract nucleotide sequences for cluster reps -----------
    //

    grep_fna = GREP_FNA (
        sequences.map { r -> tuple(r.meta, r.fna) },
        reps_stats.map { r -> r.reps_list }
    )

    gunzip_fna = UNCOMPRESSED_REPS_FNA(
        grep_fna.sequences
    )

    //
    // ----------- Extract protein sequences for cluster reps -----------
    //

    grep_faa = GREP_FAA (
        faa.map { faa_file -> tuple([id: type], faa_file) },
        reps_stats.map { r -> r.reps_proteins_list }
    )

    gunzip_faa = UNCOMPRESSED_REPS_FAA(
        grep_faa.sequences
    )

    reps = reps_stats
        .join(grep_fna.sequences.map { meta, f -> record(meta: meta, fna_gz: f) }, by: 'meta')
        .join(gunzip_fna.gunzip.map { meta, f -> record(meta: meta, fna: f) }, by: 'meta')
        .join(grep_faa.sequences.map { meta, f -> record(meta: meta, faa_gz: f) }, by: 'meta')
        .join(gunzip_faa.gunzip.map { meta, f -> record(meta: meta, faa: f) }, by: 'meta')
        .map { r ->
            record(
                meta: r.meta, reps_list: r.reps_list, reps_gff: r.reps_gff, reps_proteins_list: r.reps_proteins_list,
                fna: r.fna, faa: r.faa, fna_gz: r.fna_gz, faa_gz: r.faa_gz
            )
        }

    emit:
    reps: Channel<ClusterReps> = reps
}
