/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { SEQTK_SUBSEQ as GREP_FNA         } from '../../modules/nf-core/seqtk/subseq'
include { SEQTK_SUBSEQ as GREP_FAA         } from '../../modules/nf-core/seqtk/subseq'
include { GUNZIP as UNCOMPRESSED_REPS_FNA  } from '../../modules/nf-core/gunzip'
include { GUNZIP as UNCOMPRESSED_REPS_FAA  } from '../../modules/nf-core/gunzip'

include { BACPHLIP                         } from '../../modules/local/bacphlip'
include { CRISPRCAS_FINDER                 } from '../../modules/local/crispcasfinder'
include { EXTRACT_REPS_STATS               } from '../../modules/local/extract_reps_stats'
include { SANKEY_PLOT as SANKEY_VITAP      } from '../../modules/local/sankey_plot'
include { VITAP                            } from '../../modules/local/vitap'

include { CLUSTERING                       } from './clustering'
include { TAXONOMY_VISUALISATION           } from './taxonomy_visualisation'
include { PROTEINS_PROCESSING              } from './proteins_subwf'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PROCESS_VIRAL_SEQUENCES {

    take:
    sequences
    gff
    faa
    mapfile

    main:

    ch_versions = channel.empty()

    //
    // -------- Cluster sequences
    //
    sequences
        .filter { meta, seqs ->
            seqs && seqs.size() > 0
        }
        .set { samples_seqs }

    CLUSTERING(
       samples_seqs,
       'ani',
       false,
       false,
       0.95
    )
    ch_versions = ch_versions.mix(CLUSTERING.out.versions)

    //
    // -------- Statistics and taxonomy from GFF for viral_sequences reps
    //
    EXTRACT_REPS_STATS (
        CLUSTERING.out.clusters_tsv,
        gff,
        mapfile
    )
    ch_versions = ch_versions.mix(EXTRACT_REPS_STATS.out.versions)


    // -------- Extract sequences for cluster reps
    GREP_FNA (
        sequences,
        EXTRACT_REPS_STATS.out.reps_list.map{ id, tsv -> tsv }
    )
    ch_versions = ch_versions.mix(GREP_FNA.out.versions)

    UNCOMPRESSED_REPS_FNA(
        GREP_FNA.out.sequences
    )
    ch_versions = ch_versions.mix(UNCOMPRESSED_REPS_FNA.out.versions)


    // -------- Extract sequences for proteins cluster reps
    GREP_FAA (
        faa.map{faa_item -> [[id: 'viral_sequences'], faa_item]},
        EXTRACT_REPS_STATS.out.reps_proteins_list.map{ id, tsv -> tsv }
    )
    ch_versions = ch_versions.mix(GREP_FAA.out.versions)

    UNCOMPRESSED_REPS_FAA(
        GREP_FAA.out.sequences
    )
    ch_versions = ch_versions.mix(UNCOMPRESSED_REPS_FAA.out.versions)

    //
    // -------- Host assignment
    //
    CRISPRCAS_FINDER(
        UNCOMPRESSED_REPS_FNA.out.gunzip
    )
    ch_versions = ch_versions.mix(CRISPRCAS_FINDER.out.versions)

    //
    // -------- Lifestyle
    // predicting bacteriophage lifestyle from conserved protein domains
    //
    BACPHLIP(
        UNCOMPRESSED_REPS_FNA.out.gunzip
    )
    ch_versions = ch_versions.mix(BACPHLIP.out.versions)


    TAXONOMY_VISUALISATION(
       UNCOMPRESSED_REPS_FNA.out.gunzip,
       EXTRACT_REPS_STATS.out.reps_krona_tsv
    )

    if (params.run_vitap_taxonomy) {
        //
        // Taxonomy VITAP testing...
        //
        VITAP (
            GREP_FNA.out.sequences,
            params.vitap_db
        )
        ch_versions = ch_versions.mix(VITAP.out.versions)

        SANKEY_VITAP (
           VITAP.out.best_lineages
        )
        ch_versions = ch_versions.mix(SANKEY_VITAP.out.versions)
    }

    //
    // ----------- Proteins processing
    //
    PROTEINS_PROCESSING(
       UNCOMPRESSED_REPS_FAA.out.gunzip.join(EXTRACT_REPS_STATS.out.reps_gff)
    )
    ch_versions = ch_versions.mix(PROTEINS_PROCESSING.out.versions)

    emit:

    reps_seqs      = GREP_FNA.out.sequences  // compressed
    reps_proteins  = GREP_FAA.out.sequences  // compressed
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}
