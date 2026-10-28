nextflow.enable.types = true

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { SAMTOOLS_FAIDX as INDEX_FAA                  } from '../../../modules/nf-core/samtools/faidx'
include { SAMTOOLS_FAIDX as INDEX_FNA                  } from '../../../modules/nf-core/samtools/faidx'
include { TABIX_BGZIPTABIX as INDEX_COMPRESS_GFF       } from '../../../modules/nf-core/tabix/bgziptabix'
include { TABIX_BGZIPTABIX as BGZIP_FNA                } from '../../../modules/nf-core/tabix/bgziptabix'
include { TABIX_BGZIPTABIX as BGZIP_FAA                } from '../../../modules/nf-core/tabix/bgziptabix'
include { TABIX_BGZIPTABIX as INDEX_COMPRESS_BACPHLIP  } from '../../../modules/nf-core/tabix/bgziptabix'

include { SORT_GFF                                     } from '../../../modules/local/sort_gff'


record IndexInput {
    meta: Map
    reps_gff: Path
    fna: Path   // uncompressed
    faa: Path   // uncompressed
}

record BacphlipTable {
    meta: Map
    bacphlip: Path
}

workflow INDEX_RESULTS {

    take:
    reps: Channel<IndexInput>
    bacphlip: Channel<BacphlipTable>   // empty for plasmids

    main :

    //
    // ----------- viral reps GFF -----------
    //
    sorted_gff = SORT_GFF (
        reps.map { r -> record(meta: r.meta, gff: r.reps_gff) }
    )

    // Index and compress GFF
    INDEX_COMPRESS_GFF (
        sorted_gff.map { r -> tuple(r.meta, r.sorted_gff) }
    )

    //
    // ----------- FNA (nucleotides) -----------
    //
    bgzip_fna = BGZIP_FNA (
        reps.map { r -> tuple(r.meta, r.fna) }
    )

    INDEX_FNA (
        bgzip_fna.gz_index.map { meta, fasta, _index -> tuple(meta, fasta, []) },
        false
    )

    //
    // ----------- FAA (proteins) -----------
    //
    bgzip_faa = BGZIP_FAA (
        reps.map { r -> tuple(r.meta, r.faa) }
    )

    INDEX_FAA (
        bgzip_faa.gz_index.map { meta, fasta, _index -> tuple(meta, fasta, []) },
        false
    )

    //
    // ----------- bacphlip -----------
    //
    INDEX_COMPRESS_BACPHLIP (
        bacphlip.map { r -> tuple(r.meta, r.bacphlip) }
    )
}
