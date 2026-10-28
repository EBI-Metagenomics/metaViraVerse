nextflow.enable.types = true

include { CHOOSE_SEQUENCES                               } from '../../../modules/local/choose_sequences'
include { RENAME_CONTIGS                                 } from '../../../modules/local/rename_contigs'
include { SEPARATE_SEQUENCES as SEPARATE_VIRAL_SEQUENCES } from '../../../modules/local/separate_sequences'
include { SEPARATE_SEQUENCES as SEPARATE_PLASMIDS        } from '../../../modules/local/separate_sequences'
include { SEPARATE_SEQUENCES as SEPARATE_PROPHAGES       } from '../../../modules/local/separate_sequences'

include { BARRNAP                                        } from '../../../modules/nf-core/barrnap'
include { CSVTK_CONCAT as CONCATENATE_CHECKV             } from '../../../modules/nf-core/csvtk/concat'
include { CHECKV_ENDTOEND                                } from '../../../modules/nf-core/checkv/endtoend'
include { SEQKIT_SPLIT2 as CHUNK_FNA                     } from '../../../modules/nf-core/seqkit/split2'
include { PIGZ_COMPRESS as COMPRESS_VIRAL_SEQ_FILES      } from '../../../modules/nf-core/pigz/compress/main'
include { PIGZ_COMPRESS as COMPRESS_PLASMIDS_FILES       } from '../../../modules/nf-core/pigz/compress/main'
include { PIGZ_COMPRESS as COMPRESS_PROPHAGES_FILES      } from '../../../modules/nf-core/pigz/compress/main'


record Sample {
    meta: Map    // [id, source, biome] (MGnify) or [id, type, source, biome, ...] (third-party)
    fna: Path
    gff: Path
    faa: Path
}

record RenameParams {
    start_accession: Integer?
    end_accession: Integer?
    viral_sequence_identifier: String?
    prophage_identifier: String?
    plasmid_identifier: String?
    rename_accession: String
}

record Preprocessed {
    meta: Map                   // [id: 'combined']
    map_file: Path              // rename map
    metadata: Path              // metadata of all records
    filtered_metadata: Path
    excluded_metadata: Path
    filtered_gff: Path          // concatenated GFF built from all inputs
    filtered_faa: Path          // concatenated FAA built from all inputs
    virus_fna: Path
    virus_gff: Path
    virus_faa: Path
    prophage_fna: Path
    prophage_gff: Path
    prophage_faa: Path
    plasmid_fna: Path
    plasmid_gff: Path
    plasmid_faa: Path
    stats: Path
}

workflow PREPROCESSING {

    take:
    samples: Channel<Sample>
    third_party: Channel<Sample>
    rename_params: RenameParams

    main:

    // ----------- Assign unique identifiers to all coming sequences by each record in samplesheet
    // We will combine all sequences together and keep track of names, biomes and types in map file
    // The new names would be {prefix}{number}, ex. >seq1
    // For MGnify processing rename_accession should be MGYV
    // It will also rename ID in attributes column in GFF
    // That step is running with --combine option and it will return one renamed FASTA and GFF
    // Samples are collected once (in arrival order, which sets the accession numbering),
    // so that the per-field lists below stay aligned; collect emits an empty list when
    // there are no samples of that kind.

    def combined_meta = [id: 'combined']
    def mgnify      = samples.collect().map { ss -> ss.toList() }
    def third_party_list = third_party.collect().map { ss -> ss.toList() }

    renamed = RENAME_CONTIGS(
        combined_meta,
        mgnify.map { ss -> ss.collect { s -> s.fna }.toList() },
        mgnify.map { ss -> ss.collect { s -> s.gff }.toList() },
        mgnify.map { ss -> ss.collect { s -> s.meta.source as String }.toList() },
        mgnify.map { ss -> ss.collect { s -> s.meta.biome as String }.toList() },
        rename_params.start_accession,
        rename_params.end_accession,
        third_party_list.map { ss -> ss.collect { s -> s.fna }.toList() },
        third_party_list.map { ss -> ss.collect { s -> s.gff }.toList() },
        third_party_list.map { ss -> ss.collect { s -> s.meta.type as String }.toList() },
        third_party_list.map { ss -> ss.collect { s -> s.meta.biome as String }.toList() },
        third_party_list.map { ss -> ss.collect { s -> s.meta.source as String }.toList() },
        rename_params.viral_sequence_identifier,
        rename_params.prophage_identifier,
        rename_params.plasmid_identifier,
        rename_params.rename_accession
    )

    // Raw pool of protein sequences (original, unrenamed protein IDs) from both sources,
    // used later to pull out just the proteins that survive separation + deduplication
    def faa_raw = samples.map { r -> r.faa }
        .mix( third_party.map { r -> r.faa } )
        .collectFile( name: 'combined_raw.faa' )
        .collect()
        .map { files -> files.toList()[0] }

    //
    // ----- SEPARATE SEQUENCES INTO VIRUSES, PROPHAGES AND PLASMIDS ------
    // Splitting happens right after renaming, based on the 'definition' column of the
    // rename map (third_party_* definitions are folded into the same bucket as their
    // MGnify counterpart). FNA, GFF and FAA are all split together so downstream steps
    // (quality evaluation, choose_sequences) have a consistent triple per category.
    //
    def renamed_fna = renamed.map { r -> record(meta: r.meta, fna: r.fna_renamed) }
    def renamed_gff = renamed.map { r -> r.gff_renamed }
    def map_file    = renamed.map { r -> r.map_file }

    viruses = SEPARATE_VIRAL_SEQUENCES(
       renamed_fna,
       renamed_gff,
       faa_raw,
       map_file,
       'virus'
    )

    prophages = SEPARATE_PROPHAGES(
       renamed_fna,
       renamed_gff,
       faa_raw,
       map_file,
       'prophage'
    )

    plasmids = SEPARATE_PLASMIDS(
       renamed_fna,
       renamed_gff,
       faa_raw,
       map_file,
       'plasmid'
    )

    def ch_virus_prophage_fna = viruses.combine(prophages)
        .flatMap { v, p -> [v.fna, p.fna] }
        .collectFile( name: 'viruses_and_prophages.fna' )
        .map { fna -> tuple(combined_meta, fna) }

    //
    // ----------- Evaluate quality and rRNA content for viruses + prophages only
    // (plasmids are not checked for viral quality/rRNA contamination)
    //
    def quality: Value<Path?> = channel.value(null)
    if ( ! params.skip_checkv ) {
        chunks = CHUNK_FNA (
            ch_virus_prophage_fna,
            [],                                        // length: (disabled) max number of nucleotides per chunk
            params.nucleotide_fasta_chunksize,         // size: max number of sequences per chunk
        )
        // a single chunk is emitted as a path, several as a list: transpose handles both
        def ch_fna_chunks = chunks.chunked_output.transpose()

        checkv = CHECKV_ENDTOEND (
            ch_fna_chunks,
            params.checkv_db
        )

        checkv_quality = CONCATENATE_CHECKV (
            checkv.quality_summary.groupTuple(),
            'tsv',
            'tsv'
        )
        quality = checkv_quality.csv
            .map { _meta, table -> table }
            .collect()
            .map { tables ->
                if( tables.isEmpty() )
                    error "CheckV produced no quality summary -- refusing to choose sequences without quality filtering"
                tables.toList()[0]
            }
    }

    //
    // ----- Detect rRNA sequences ------
    //
    def rna_gff: Value<Path?> = channel.value(null)
    if ( ! params.skip_rrna_detection ) {
        barrnap = BARRNAP(
          ch_virus_prophage_fna.map { meta, fasta -> tuple(meta, fasta, "bac") }
        )

        rna_gff = barrnap.gff
            .map { _meta, gff -> gff }
            .collect()
            .map { gffs ->
                if( gffs.isEmpty() )
                    error "barrnap produced no GFF -- refusing to choose sequences without rRNA filtering"
                gffs.toList()[0]
            }
    }

    //
    // ----------- Deduplicate across categories, leave unique and save metadata
    // Compares hashsums across viruses/prophages/plasmids together (in case the same
    // sequence was classified differently by different predictions), prioritising
    // metagenome over genome and MGnify over third-party.
    // Also filters out non-determined-quality viruses/prophages, and drops the
    // corresponding GFF/FAA records for any sequence that gets deduplicated away.
    //
    chosen = CHOOSE_SEQUENCES (
        renamed.map { r -> record(meta: r.meta, map_file: r.map_file) },
        viruses,
        prophages,
        plasmids,
        quality,
        rna_gff
    )

    COMPRESS_VIRAL_SEQ_FILES (
        chosen
            .flatMap { r -> [r.virus_fna, r.virus_gff, r.virus_faa].collect { f -> tuple(r.meta, f) } }
    )

    COMPRESS_PLASMIDS_FILES (
        chosen
            .flatMap { r -> [r.plasmid_fna, r.plasmid_gff, r.plasmid_faa].collect { f -> tuple(r.meta, f) } }
    )

    COMPRESS_PROPHAGES_FILES (
        chosen
            .flatMap { r -> [r.prophage_fna, r.prophage_gff, r.prophage_faa].collect { f -> tuple(r.meta, f) } }
    )

    preprocessed = chosen
        .map { r ->
            record(
                meta: r.meta,
                metadata: r.metadata, filtered_metadata: r.filtered_metadata, excluded_metadata: r.excluded_metadata,
                filtered_gff: r.filtered_gff, filtered_faa: r.filtered_faa,
                virus_fna: r.virus_fna, virus_gff: r.virus_gff, virus_faa: r.virus_faa,
                prophage_fna: r.prophage_fna, prophage_gff: r.prophage_gff, prophage_faa: r.prophage_faa,
                plasmid_fna: r.plasmid_fna, plasmid_gff: r.plasmid_gff, plasmid_faa: r.plasmid_faa,
                stats: r.stats
            )
        }
        .combine(map_file: map_file)

    emit:
    preprocessed: Value<Preprocessed> = preprocessed
}
