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

workflow EXTRACT_CLUSTER_FILES {

    take:
    clusters_tsv
    gff
    mapfile
    sequences
    faa
    type

    main:
    ch_versions = channel.empty()

    //
    // ----------- Statistics and taxonomy from GFF for reps -----------
    //

    EXTRACT_REPS_STATS (
        clusters_tsv,
        gff,
        mapfile
    )
    ch_versions = ch_versions.mix(EXTRACT_REPS_STATS.out.versions)

    //
    // ----------- Extract nucleotide sequences for cluster reps -----------
    //

    GREP_FNA (
        sequences,
        EXTRACT_REPS_STATS.out.reps_list.map{ id, tsv -> tsv }
    )
    ch_versions = ch_versions.mix(GREP_FNA.out.versions)

    UNCOMPRESSED_REPS_FNA(
        GREP_FNA.out.sequences
    )

    //
    // ----------- Extract protein sequences for cluster reps -----------
    //

    GREP_FAA (
        faa.map{faa_item -> [[id: type], faa_item]},
        EXTRACT_REPS_STATS.out.reps_proteins_list.map{ id, tsv -> tsv }
    )
    ch_versions = ch_versions.mix(GREP_FAA.out.versions)

    UNCOMPRESSED_REPS_FAA(
        GREP_FAA.out.sequences
    )

    emit:
    reps_fna_uncompressed      = UNCOMPRESSED_REPS_FNA.out.gunzip
    reps_faa_uncompressed      = UNCOMPRESSED_REPS_FAA.out.gunzip
    reps_fna_compressed        = GREP_FNA.out.sequences
    reps_faa_compressed        = GREP_FAA.out.sequences
    reps_gff                   = EXTRACT_REPS_STATS.out.reps_gff
    reps_list                  = EXTRACT_REPS_STATS.out.reps_list
    versions                   = ch_versions
}
