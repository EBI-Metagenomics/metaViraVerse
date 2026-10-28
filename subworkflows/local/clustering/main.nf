nextflow.enable.types = true

include { ANICALC               } from '../../../modules/local/checkv/anicalc'
include { ANICLUST              } from '../../../modules/local/checkv/aniclust'

include { BLAST_MAKEBLASTDB     } from '../../../modules/nf-core/blast/makeblastdb'
include { BLAST_BLASTN          } from '../../../modules/nf-core/blast/blastn'
include { VCLUST_ALIGN          } from '../../../modules/nf-core/vclust/align'
include { VCLUST_CLUSTER        } from '../../../modules/nf-core/vclust/cluster'
include { VCLUST_PREFILTER      } from '../../../modules/nf-core/vclust/prefilter'


record ClusteringInput {
    meta: Map
    fna: Path
}

record Clusters {
    meta: Map
    clusters_tsv: Path
}

workflow CLUSTERING {

    take:
    sequences: Channel<ClusteringInput>
    vclust_metric: String
    // thresholds: 0 = unset (the untyped VCLUST_CLUSTER module skips falsy values and rejects null)
    vclust_tani_threshold: Float
    vclust_gani_threshold: Float
    vclust_ani_threshold: Float
    cluster_vclust: Boolean   // vclust (true) or blastn + CheckV ANI clustering (false)

    main:

    def ch_clusters: Channel<Clusters> = channel.empty()
    def ch_fna = sequences.map { r -> tuple(r.meta, r.fna) }

    if ( cluster_vclust ) {
        prefilter = VCLUST_PREFILTER (
            ch_fna
        )

        align = VCLUST_ALIGN (
            ch_fna,
            prefilter.txt,
            false
        )

        vclust = VCLUST_CLUSTER (
            align.tsv,
            align.ids,
            vclust_metric,
            vclust_tani_threshold,
            vclust_gani_threshold,
            vclust_ani_threshold
        )

        ch_clusters = vclust.clusters.map { meta, tsv -> record(meta: meta, clusters_tsv: tsv) }

    } else {
        // Creation of a blast+ database
        blastdb = BLAST_MAKEBLASTDB(
            ch_fna,
            'nucl'
        )

        // Using megablast from blast+ package to perform all-vs-all blastn of sequences
        blastn = BLAST_BLASTN(
            ch_fna,
            blastdb.db,
            [],
            '',
            ''
        )

        // Calculate pairwise ANI by combining local alignments between sequence pairs
        ani = ANICALC(
           blastn.txt.map { meta, tsv -> record(meta: meta, tsv: tsv) }
        )

        // UCLUST-like clustering
        // simple greedy, centroid-based clustering
        // sometimes called leader clustering or single-pass centroid clustering
        ch_clusters = ANICLUST(
           ani.join(sequences, by: 'meta')
               .map { r -> record(meta: r.meta, ani_tsv: r.calculated_ani_tsv, sequences: r.fna) }
        )
    }

    clusters = ch_clusters

    emit:
    clusters: Channel<Clusters> = clusters

}
