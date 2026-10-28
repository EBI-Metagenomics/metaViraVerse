nextflow.enable.types = true

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FIND_CONCATENATE as CONCATENATE_IPHOP_GENOME     } from '../../../modules/nf-core/find/concatenate'
include { FIND_CONCATENATE as CONCATENATE_IPHOP_GENUS      } from '../../../modules/nf-core/find/concatenate'
include { IPHOP_PREDICT                                    } from '../../../modules/nf-core/iphop/predict/main'
include { SEQKIT_SPLIT2 as CHUNK_FNA_IPHOP                 } from '../../../modules/nf-core/seqkit/split2'

include { COLLECT_HOST_INFO                                } from '../../../modules/local/spacepharer/collect_crispr_host_info'
include { SPACEPHARER_CREATEDB as SPACEPHARER_CREATEDB     } from '../../../modules/local/spacepharer/createdb'
include { SPACEPHARER_CREATEDB as SPACEPHARER_CREATEDB_REV } from '../../../modules/local/spacepharer/createdb'
include { SPACEPHARER_EASYPREDICT                          } from '../../../modules/local/spacepharer/easy_predict'
include { CHANGE_SPACE_TO_UNDERSCORE                       } from '../../../modules/local/change_space_to_underscore'
include { SANKEY_PLOT as SANKEY_IPHOP_GENUS_HOST           } from '../../../modules/local/sankey_plot/main'
include { SANKEY_PLOT as SANKEY_IPHOP_GENOME_HOST          } from '../../../modules/local/sankey_plot/main'
include { SANKEY_PLOT as SANKEY_IPHOP_GENUS_HOST_REDUCED   } from '../../../modules/local/sankey_plot/main'
include { SANKEY_PLOT as SANKEY_IPHOP_GENOME_HOST_REDUCED  } from '../../../modules/local/sankey_plot/main'
include { SANKEY_PLOT as SANKEY_SPACEPHARER_HOST           } from '../../../modules/local/sankey_plot/main'
include { SANKEY_PLOT as SANKEY_SPACEPHARER_HOST_REDUCED   } from '../../../modules/local/sankey_plot/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

record HostInput {
    meta: Map
    fna: Path   // representative sequences, uncompressed
}

record IphopGenome {
    meta: Map
    iphop_genome: Path
}

record IphopGenus {
    meta: Map
    iphop_genus: Path
}

record SpacepharerHost {
    meta: Map
    spacepharer: Path
}

record HostResult {
    meta: Map
    iphop_genome: Path?   // null when iPHoP is skipped
    iphop_genus: Path?    // null when iPHoP is skipped
    spacepharer: Path?    // null without --predict_host_from_custom_spacers
}

workflow HOST_DETECTION {

    take:
    fna: Channel<HostInput>
    metadata: String?       // --custom_genomes_metadata
    mapfile: Value<Path>

    main:

    def iphop_host_genome: Channel<IphopGenome> = channel.empty()
    def iphop_host_genus: Channel<IphopGenus> = channel.empty()
    def spacepharer_host: Channel<SpacepharerHost> = channel.empty()

    if ( !params.skip_iphop) {
        chunks = CHUNK_FNA_IPHOP (
            fna.map { r -> tuple(r.meta, r.fna) },
            [],                                        // length: (disabled) max number of nucleotides per chunk
            params.nucleotide_fasta_chunksize_iphop,   // size: max number of sequences per chunk
        )
        // a single chunk is emitted as a path, several as a list: transpose handles both
        def ch_fna_chunks = chunks.chunked_output.transpose()

        iphop = IPHOP_PREDICT (
            ch_fna_chunks,
            params.iphop_db
        )

        iphop_genome = CONCATENATE_IPHOP_GENOME (
            iphop.iphop_genome
                .map { meta, files -> tuple(meta + [id: "${meta.id}_genome"], files) }
                .groupTuple(),
            1
        )

        iphop_genus = CONCATENATE_IPHOP_GENUS (
            iphop.iphop_genus
                .map { meta, files -> tuple(meta + [id: "${meta.id}_genus"], files) }
                .groupTuple(),
            1
        )

        // Plot Sankey for hosts reduced to genus and top 20 hits for website
        SANKEY_IPHOP_GENUS_HOST_REDUCED (
            iphop_genus.file_out.map { meta, table -> record(meta: meta, table: table) }
        )
        SANKEY_IPHOP_GENOME_HOST_REDUCED (
            iphop_genome.file_out.map { meta, table -> record(meta: meta, table: table) }
        )
        // Plot Sankey for hosts
        SANKEY_IPHOP_GENUS_HOST (
            iphop_genus.file_out.map { meta, table -> record(meta: meta, table: table) }
        )
        SANKEY_IPHOP_GENOME_HOST (
            iphop_genome.file_out.map { meta, table -> record(meta: meta, table: table) }
        )

        // the tables carry a renamed meta (<id>_genome / <id>_genus): attach the input's meta
        iphop_host_genome = fna
            .combine(iphop_genome.file_out.map { _meta, table -> record(iphop_genome: table) })
            .map { r, host -> record(meta: r.meta, iphop_genome: host.iphop_genome) }
        iphop_host_genus = fna
            .combine(iphop_genus.file_out.map { _meta, table -> record(iphop_genus: table) })
            .map { r, host -> record(meta: r.meta, iphop_genus: host.iphop_genus) }
    }

    if (params.predict_host_from_custom_spacers) {
        sanitised = CHANGE_SPACE_TO_UNDERSCORE(
            fna.map { r -> record(meta: r.meta, fasta: r.fna) }
        )
        sanitised_fasta = sanitised.map { r -> record(meta: r.meta, fasta: r.sanitised_fasta) }

        target_db = SPACEPHARER_CREATEDB (
            sanitised_fasta,
            "targetSetDB"
        )

        target_db_rev = SPACEPHARER_CREATEDB_REV (
            sanitised_fasta,
            "targetSetDB_rev"
        )

        predictions = SPACEPHARER_EASYPREDICT(
            target_db.map { r -> record(meta: r.meta, target_db: r.db) }
                .join(target_db_rev.map { r -> record(meta: r.meta, target_db_rev: r.db) }, by: 'meta')
                .map { r ->
                    record(
                        meta: [id: 'spacers'],
                        crispr_fasta: file(params.custom_spacers_fasta),
                        target_db: r.target_db,
                        target_db_rev: r.target_db_rev
                    )
                },
            "targetSetDB"
        )

        host_info = COLLECT_HOST_INFO (
           predictions,
           file(params.custom_spacers_metadata),
           mapfile,
           file(metadata)
        )

        SANKEY_SPACEPHARER_HOST (
            host_info.map { r -> record(meta: r.meta, table: r.tsv) }
        )

        SANKEY_SPACEPHARER_HOST_REDUCED (
            host_info.map { r -> record(meta: r.meta, table: r.tsv) }
        )

        // predictions carry the spacers' meta: attach the input's meta
        spacepharer_host = fna
            .combine(host_info)
            .map { r, host -> record(meta: r.meta, spacepharer: host.tsv) }
    }

    // Steps above are optional: their fields are null when skipped
    hosts = fna
        .map { r -> record(meta: r.meta) }
        .join(iphop_host_genome, by: 'meta', remainder: true)
        .join(iphop_host_genus, by: 'meta', remainder: true)
        .join(spacepharer_host, by: 'meta', remainder: true)
        .map { r -> record(meta: r.meta, iphop_genome: r.iphop_genome, iphop_genus: r.iphop_genus, spacepharer: r.spacepharer) }

    emit:
    hosts: Channel<HostResult> = hosts
}
