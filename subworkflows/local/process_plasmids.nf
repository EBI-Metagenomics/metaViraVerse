/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { EXTRACT_REPS_STATS                      } from '../../modules/local/extract_reps_stats'

include { SEQTK_SUBSEQ as GREP_FAA                } from '../../modules/nf-core/seqtk/subseq'
include { SEQTK_SUBSEQ as GREP_FNA                } from '../../modules/nf-core/seqtk/subseq'
include { TABIX_BGZIPTABIX as INDEX_COMPRESS_GFF  } from '../../modules/nf-core/tabix/bgziptabix'
include { TABIX_TABIX as INDEX_FAA                } from '../../modules/nf-core/tabix/tabix'
include { TABIX_TABIX as INDEX_FNA                } from '../../modules/nf-core/tabix/tabix'

include { CLUSTERING                              } from './clustering'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PROCESS_PLASMIDS {

    take:
    sequences
    combined_faa
    combined_gff
    mapfile

    main:

    ch_versions = channel.empty()

    //
    // Cluster sequences
    //
    sequences
        .filter { _meta, seqs ->
            seqs && seqs.size() > 0
        }
        .set { samples_seqs }

    CLUSTERING(
       samples_seqs,
       'gani',
       false,
       0.35,
       false
    )
    ch_versions = ch_versions.mix(CLUSTERING.out.versions)

    //
    // -------- Statistics and taxonomy from GFF for plasmid reps
    //
    EXTRACT_REPS_STATS (
        CLUSTERING.out.clusters_tsv,
        combined_gff,
        mapfile
    )
    ch_versions = ch_versions.mix(EXTRACT_REPS_STATS.out.versions)

    // Index and compress GFF
    INDEX_COMPRESS_GFF (
        EXTRACT_REPS_STATS.out.reps_gff
    )

    // Extract sequences for cluster reps
    GREP_FNA (
        sequences,
        EXTRACT_REPS_STATS.out.reps_list.map{ id, tsv -> tsv }
    )
    ch_versions = ch_versions.mix(GREP_FNA.out.versions)

    INDEX_FNA (
        GREP_FNA.out.sequences
    )

    // Extract proteins for cluster reps
    GREP_FAA (
        combined_faa.map{faa_item -> [[id: 'plasmids'], faa_item]},
        EXTRACT_REPS_STATS.out.reps_proteins_list.map{ id, tsv -> tsv }
    )
    ch_versions = ch_versions.mix(GREP_FAA.out.versions)

    INDEX_FAA (
        GREP_FAA.out.sequences
    )

    emit:

    reps_tsv       = EXTRACT_REPS_STATS.out.reps_list
    reps_seqs      = GREP_FNA.out.sequences      // compressed
    reps_proteins  = GREP_FAA.out.sequences      // compressed
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}
