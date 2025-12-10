include { ANICALC               } from '../../modules/local/checkv/anicalc'
include { ANICLUST              } from '../../modules/local/checkv/aniclust'

include { BLAST_MAKEBLASTDB     } from '../../modules/nf-core/blast/makeblastdb'
include { BLAST_BLASTN          } from '../../modules/nf-core/blast/blastn'

workflow CLUSTERING {

    take:
    sequences
    min_ani
    min_coverage

    main:

    ch_versions = Channel.empty()

    // Creation of a blast+ database
    BLAST_MAKEBLASTDB(
        sequences
    )
    ch_versions = ch_versions.mix(BLAST_MAKEBLASTDB.out.versions)


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
    ANICLUST(
       ANICALC.out.calculated_ani_tsv.join(sequences),
       min_ani,
       min_coverage
    )
    ch_versions = ch_versions.mix(ANICLUST.out.versions)

    emit:
    clusters_tsv   = ANICLUST.out.clusters_tsv
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}
