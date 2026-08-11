include { ANICALC               } from '../../../modules/local/checkv/anicalc'
include { ANICLUST              } from '../../../modules/local/checkv/aniclust'

include { BLAST_MAKEBLASTDB     } from '../../../modules/nf-core/blast/makeblastdb'
include { BLAST_BLASTN          } from '../../../modules/nf-core/blast/blastn'
include { VCLUST_ALIGN          } from '../../../modules/nf-core/vclust/align'
include { VCLUST_CLUSTER        } from '../../../modules/nf-core/vclust/cluster'
include { VCLUST_PREFILTER      } from '../../../modules/nf-core/vclust/prefilter'


workflow CLUSTERING {

    take:
    sequences
    vclust_metric
    vclust_tani_threshold
    vclust_gani_threshold
    vclust_ani_threshold

    main:

    ch_versions = channel.empty()

    if ( params.cluster_vclust ) {
        VCLUST_PREFILTER (
            sequences
        )
        ch_versions = ch_versions.mix(VCLUST_PREFILTER.out.versions)

        VCLUST_ALIGN (
            sequences,
            VCLUST_PREFILTER.out.txt,
            false
        )
        ch_versions = ch_versions.mix(VCLUST_ALIGN.out.versions)

        VCLUST_CLUSTER (
            VCLUST_ALIGN.out.tsv,
            VCLUST_ALIGN.out.ids,
            vclust_metric,
            vclust_tani_threshold,
            vclust_gani_threshold,
            vclust_ani_threshold
        )
        ch_versions = ch_versions.mix(VCLUST_CLUSTER.out.versions)

        clusters_tsv = VCLUST_CLUSTER.out.clusters

    } else {
        // Creation of a blast+ database
        BLAST_MAKEBLASTDB(
            sequences,
            'nucl'
        )

        // Using megablast from blast+ package to perform all-vs-all blastn of sequences
        BLAST_BLASTN(
            sequences,
            BLAST_MAKEBLASTDB.out.db,
            [],
            '',
            ''
        )
        ch_versions = ch_versions.mix(BLAST_BLASTN.out.versions)

        // Calculate pairwise ANI by combining local alignments between sequence pairs
        ANICALC(
           BLAST_BLASTN.out.txt
        )
        ch_versions = ch_versions.mix(ANICALC.out.versions)

        // UCLUST-like clustering
        // simple greedy, centroid-based clustering
        // sometimes called leader clustering or single-pass centroid clustering
        ANICLUST(
           ANICALC.out.calculated_ani_tsv.join(sequences)
        )
        ch_versions = ch_versions.mix(ANICLUST.out.versions)

        clusters_tsv = ANICLUST.out.clusters_tsv
    }

    emit:
    clusters_tsv   = clusters_tsv
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}
