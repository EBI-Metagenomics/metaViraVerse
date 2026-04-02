/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { GUNZIP as UNCOMPRESSED_REPS_FNA              } from '../../modules/nf-core/gunzip'
include { GUNZIP as UNCOMPRESSED_REPS_FAA              } from '../../modules/nf-core/gunzip'
include { SEQTK_SUBSEQ as GREP_FNA                     } from '../../modules/nf-core/seqtk/subseq'
include { SEQTK_SUBSEQ as GREP_FAA                     } from '../../modules/nf-core/seqtk/subseq'
include { TABIX_BGZIPTABIX as INDEX_COMPRESS_BACPHLIP  } from '../../modules/nf-core/tabix/bgziptabix'
include { TABIX_BGZIPTABIX as INDEX_COMPRESS_GFF       } from '../../modules/nf-core/tabix/bgziptabix'
include { TABIX_BGZIPTABIX as BGZIP_FNA                } from '../../modules/nf-core/tabix/bgziptabix'
include { TABIX_BGZIPTABIX as BGZIP_FAA                } from '../../modules/nf-core/tabix/bgziptabix'
include { SAMTOOLS_FAIDX as INDEX_FAA                  } from '../../modules/nf-core/samtools/faidx'
include { SAMTOOLS_FAIDX as INDEX_FNA                  } from '../../modules/nf-core/samtools/faidx'

include { BACPHLIP                                     } from '../../modules/local/bacphlip'
include { BUILD_FINAL_GFF                              } from '../../modules/local/build_final_gff'
include { CRISPRCAS_FINDER                             } from '../../modules/local/crispcasfinder'
include { EXTRACT_REPS_STATS                           } from '../../modules/local/extract_reps_stats'
include { GENERATE_TAXONOMY_TABLE as TAX_VIPHOGS       } from '../../modules/local/generate_taxonomy_table'
include { GENERATE_TAXONOMY_TABLE as TAX_VITAP         } from '../../modules/local/generate_taxonomy_table'
include { VITAP                                        } from '../../modules/local/vitap'
include { SORT_GFF                                     } from '../../modules/local/sort_gff'

include { CLUSTERING                                   } from './clustering'
include { TAXONOMY_VISUALISATION as VIS_VIPHOGS        } from './taxonomy_visualisation'
include { TAXONOMY_VISUALISATION as VIS_VITAP          } from './taxonomy_visualisation'
include { PROTEINS_PROCESSING                          } from './proteins_subwf'


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
    combined_metadata
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
    // -------- Statistics and taxonomy from GFF for viruses reps
    //
    EXTRACT_REPS_STATS (
        CLUSTERING.out.clusters_tsv,
        gff,
        mapfile
    )
    ch_versions = ch_versions.mix(EXTRACT_REPS_STATS.out.versions)

    // --- reps GFF
    SORT_GFF (
        EXTRACT_REPS_STATS.out.reps_gff
    )
    // Index and compress GFF
    INDEX_COMPRESS_GFF (
        SORT_GFF.out.sorted_gff
    )

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

    BGZIP_FNA (UNCOMPRESSED_REPS_FNA.out.gunzip)

    INDEX_FNA (
        BGZIP_FNA.out.gz_index.map{ meta, fasta, index -> [meta, fasta, []] },
        false
    )

    // -------- Extract sequences for proteins cluster reps
    GREP_FAA (
        faa.map{faa_item -> [[id: 'viruses'], faa_item]},
        EXTRACT_REPS_STATS.out.reps_proteins_list.map{ id, tsv -> tsv }
    )
    ch_versions = ch_versions.mix(GREP_FAA.out.versions)

    UNCOMPRESSED_REPS_FAA(
        GREP_FAA.out.sequences
    )
    ch_versions = ch_versions.mix(UNCOMPRESSED_REPS_FAA.out.versions)

    BGZIP_FAA (UNCOMPRESSED_REPS_FAA.out.gunzip)

    INDEX_FAA (
        BGZIP_FAA.out.gz_index.map{ meta, fasta, index -> [meta, fasta, []] },
        false
    )

    //
    // -------- Host assignment
    //
    if (params.run_crisprcasfinder) {
        CRISPRCAS_FINDER(
            UNCOMPRESSED_REPS_FNA.out.gunzip
        )
        ch_versions = ch_versions.mix(CRISPRCAS_FINDER.out.versions)
    }

    //
    // -------- Lifestyle
    // predicting bacteriophage lifestyle from conserved protein domains
    //
    BACPHLIP(
        UNCOMPRESSED_REPS_FNA.out.gunzip
    )
    ch_versions = ch_versions.mix(BACPHLIP.out.versions)

    INDEX_COMPRESS_BACPHLIP (
        BACPHLIP.out.bacphlip_table
    )

    //
    // Taxonomy ViPhOGs
    //
    TAX_VIPHOGS (
        EXTRACT_REPS_STATS.out.reps_stats_tsv
        .map { meta, table ->
            def new_meta = meta.clone()
            new_meta.tool = 'viphogs'
            tuple(new_meta, table)
        },
        combined_metadata
    )

    VIS_VIPHOGS(
       TAX_VIPHOGS.out.taxonomy_and_metadata
    )

    //
    // Taxonomy VITAP
    //
    VITAP (
        GREP_FNA.out.sequences,
        params.vitap_db
    )
    ch_versions = ch_versions.mix(VITAP.out.versions)

    TAX_VITAP (
        VITAP.out.best_lineages.map { meta, table ->
            def new_meta = meta.clone()
            new_meta.tool = 'vitap'
            tuple(new_meta, table)
        },
        combined_metadata
    )

    VIS_VITAP (
       TAX_VITAP.out.taxonomy_and_metadata
    )

    //
    // ----------- Proteins processing
    //
    PROTEINS_PROCESSING(
       UNCOMPRESSED_REPS_FAA.out.gunzip.join(EXTRACT_REPS_STATS.out.reps_gff)
    )
    ch_versions = ch_versions.mix(PROTEINS_PROCESSING.out.versions)

    //
    // --- Build final aggregated GFF
    //
    BUILD_FINAL_GFF (
        EXTRACT_REPS_STATS.out.reps_gff,
        BACPHLIP.out.bacphlip_table,
        PROTEINS_PROCESSING.out.hmmer_tables,
        PROTEINS_PROCESSING.out.amr_gff
    )
    ch_versions = ch_versions.mix(BUILD_FINAL_GFF.out.versions)

    emit:

    reps_tsv       = EXTRACT_REPS_STATS.out.reps_list
    reps_seqs      = GREP_FNA.out.sequences  // compressed
    reps_proteins  = GREP_FAA.out.sequences  // compressed
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}
