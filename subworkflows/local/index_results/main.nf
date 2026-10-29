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

record IndexedReps {
    meta: Map
    fna_gz: Path
    fna_fai: Path?
    fna_gzi: Path?
    faa_gz: Path
    faa_fai: Path?
    faa_gzi: Path?
    gff_gz: Path
    gff_tbi: Path            // .tbi or .csi
    bacphlip_gz: Path?       // viruses only
    bacphlip_tbi: Path?
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
    gff_index = INDEX_COMPRESS_GFF (
        sorted_gff.map { r -> tuple(r.meta, r.sorted_gff) }
    )

    //
    // ----------- FNA (nucleotides) -----------
    //
    bgzip_fna = BGZIP_FNA (
        reps.map { r -> tuple(r.meta, r.fna) }
    )

    fna_index = INDEX_FNA (
        bgzip_fna.gz_index.map { meta, fasta, _index -> tuple(meta, fasta, []) },
        false
    )

    //
    // ----------- FAA (proteins) -----------
    //
    bgzip_faa = BGZIP_FAA (
        reps.map { r -> tuple(r.meta, r.faa) }
    )

    faa_index = INDEX_FAA (
        bgzip_faa.gz_index.map { meta, fasta, _index -> tuple(meta, fasta, []) },
        false
    )

    //
    // ----------- bacphlip -----------
    //
    bacphlip_index = INDEX_COMPRESS_BACPHLIP (
        bacphlip.map { r -> tuple(r.meta, r.bacphlip) }
    )

    // samtools faidx indices are optional outputs: join with remainder, missing ones are null
    indexed = gff_index.gz_index.map { meta, gz, idx -> record(meta: meta, gff_gz: gz, gff_tbi: idx) }
        .join(bgzip_fna.gz_index.map { meta, gz, _idx -> record(meta: meta, fna_gz: gz) }, by: 'meta')
        .join(bgzip_faa.gz_index.map { meta, gz, _idx -> record(meta: meta, faa_gz: gz) }, by: 'meta')
        .join(fna_index.fai.map { meta, f -> record(meta: meta, fna_fai: f) }, by: 'meta', remainder: true)
        .join(fna_index.gzi.map { meta, f -> record(meta: meta, fna_gzi: f) }, by: 'meta', remainder: true)
        .join(faa_index.fai.map { meta, f -> record(meta: meta, faa_fai: f) }, by: 'meta', remainder: true)
        .join(faa_index.gzi.map { meta, f -> record(meta: meta, faa_gzi: f) }, by: 'meta', remainder: true)
        .join(bacphlip_index.gz_index.map { meta, gz, idx -> record(meta: meta, bacphlip_gz: gz, bacphlip_tbi: idx) }, by: 'meta', remainder: true)
        .map { r ->
            record(
                meta: r.meta,
                fna_gz: r.fna_gz, fna_fai: r.fna_fai, fna_gzi: r.fna_gzi,
                faa_gz: r.faa_gz, faa_fai: r.faa_fai, faa_gzi: r.faa_gzi,
                gff_gz: r.gff_gz, gff_tbi: r.gff_tbi,
                bacphlip_gz: r.bacphlip_gz, bacphlip_tbi: r.bacphlip_tbi
            )
        }

    emit:
    indexed: Channel<IndexedReps> = indexed
}
