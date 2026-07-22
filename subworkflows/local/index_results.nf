/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { SAMTOOLS_FAIDX as INDEX_FAA                  } from '../../modules/nf-core/samtools/faidx'
include { SAMTOOLS_FAIDX as INDEX_FNA                  } from '../../modules/nf-core/samtools/faidx'
include { TABIX_BGZIPTABIX as INDEX_COMPRESS_GFF       } from '../../modules/nf-core/tabix/bgziptabix'
include { TABIX_BGZIPTABIX as BGZIP_FNA                } from '../../modules/nf-core/tabix/bgziptabix'
include { TABIX_BGZIPTABIX as BGZIP_FAA                } from '../../modules/nf-core/tabix/bgziptabix'
include { TABIX_BGZIPTABIX as INDEX_COMPRESS_BACPHLIP  } from '../../modules/nf-core/tabix/bgziptabix'

include { SORT_GFF                                     } from '../../modules/local/sort_gff'


workflow INDEX_RESULTS {

    take:
    reps_gff
    fna_uncompressed
    faa_uncompressed
    bacphlip_concatenated

    main :

    //
    // ----------- viral reps GFF -----------
    //
    SORT_GFF (
        reps_gff
    )

    // Index and compress GFF
    INDEX_COMPRESS_GFF (
        SORT_GFF.out.sorted_gff
    )

    //
    // ----------- FNA (nucleotides) -----------
    //
    BGZIP_FNA (fna_uncompressed)

    INDEX_FNA (
        BGZIP_FNA.out.gz_index.map{ meta, fasta, index -> [meta, fasta, []] },
        false
    )

    //
    // ----------- FAA (proteins) -----------
    //
    BGZIP_FAA (faa_uncompressed)

    INDEX_FAA (
        BGZIP_FAA.out.gz_index.map{ meta, fasta, index -> [meta, fasta, []] },
        false
    )

    //
    // ----------- bacphlip -----------
    //
    if ( bacphlip_concatenated ) {
        INDEX_COMPRESS_BACPHLIP (
            bacphlip_concatenated
        )
    }
}
