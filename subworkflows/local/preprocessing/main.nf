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


workflow PREPROCESSING {

    take:
    input   // [meta, gff, fna, faa]
    ch_third_party_data  // [meta[id, biome, type, source], fasta, gff, faa]

    main:

    // ----------- Assign unique identifiers to all coming sequences by each record in samplesheet
    // We will combine all sequences together and keep track of names, biomes and types in map file
    // The new names would be {prefix}{number}, ex. >seq1
    // For MGnify processing rename_accession should be MGYV
    // It will also rename ID in attributes column in GFF
    // That step is running with --combine option and it will return one renamed FASTA and GFF

    ch_fna       = input.map { meta, fna, gff, faa -> fna }.collect().ifEmpty([])
    ch_gff       = input.map { meta, fna, gff, faa -> gff }.collect().ifEmpty([])
    ch_sources   = input.map { meta, fna, gff, faa -> tuple([meta.source]) }.collect().ifEmpty([])
    ch_biomes    = input.map { meta, fna, gff, faa -> tuple([meta.biome]) }.collect().ifEmpty([])

    ch_fna_tp    = ch_third_party_data.map { meta, fna, gff, faa -> fna }.collect().ifEmpty([])
    ch_gff_tp    = ch_third_party_data.map { meta, fna, gff, faa -> gff }.collect().ifEmpty([])
    ch_types_tp  = ch_third_party_data.map { meta, fna, gff, faa -> tuple([meta.type]) }.collect().ifEmpty([])
    ch_biomes_tp = ch_third_party_data.map { meta, fna, gff, faa -> tuple([meta.biome]) }.collect().ifEmpty([])
    ch_source_tp = ch_third_party_data.map { meta, fna, gff, faa -> tuple([meta.source]) }.collect().ifEmpty([])

    // Raw pool of protein sequences (original, unrenamed protein IDs) from both sources,
    // used later to pull out just the proteins that survive separation + deduplication
    def combined_meta = [id: 'combined']

    ch_faa_raw = input.map { meta, fna, gff, faa -> faa }
        .mix( ch_third_party_data.map { meta, fna, gff, faa -> faa } )
        .collectFile( name: 'combined_raw.faa' )
        .map { faa -> [combined_meta, faa] }

    RENAME_CONTIGS(
        combined_meta,
        ch_fna,
        ch_gff,
        ch_sources,
        ch_biomes,
        params.start_accession,
        params.end_accession,
        ch_fna_tp,
        ch_gff_tp,
        ch_types_tp,
        ch_biomes_tp,
        ch_source_tp,
        params.viral_sequence_identifier,
        params.prophage_identifier,
        params.plasmid_identifier,
        params.rename_accession
    )
    ch_fna_renamed = RENAME_CONTIGS.out.map { r -> [r.meta, r.fna_renamed] }
    ch_gff_renamed = RENAME_CONTIGS.out.map { r -> r.gff_renamed }
    ch_map_file    = RENAME_CONTIGS.out.map { r -> [r.meta, r.map_file] }

    //
    // ----- SEPARATE SEQUENCES INTO VIRUSES, PROPHAGES AND PLASMIDS ------
    // Splitting happens right after renaming, based on the 'definition' column of the
    // rename map (third_party_* definitions are folded into the same bucket as their
    // MGnify counterpart). FNA, GFF and FAA are all split together so downstream steps
    // (quality evaluation, choose_sequences) have a consistent triple per category.
    //

    SEPARATE_VIRAL_SEQUENCES(
       ch_fna_renamed.map { meta, fna -> record(meta: meta, fna: fna) },
       ch_gff_renamed,
       ch_faa_raw.map { _meta, faa -> faa },
       ch_map_file.map { _meta, map_file -> map_file },
       'virus'
    )

    SEPARATE_PROPHAGES(
       ch_fna_renamed.map { meta, fna -> record(meta: meta, fna: fna) },
       ch_gff_renamed,
       ch_faa_raw.map { _meta, faa -> faa },
       ch_map_file.map { _meta, map_file -> map_file },
       'prophage'
    )

    SEPARATE_PLASMIDS(
       ch_fna_renamed.map { meta, fna -> record(meta: meta, fna: fna) },
       ch_gff_renamed,
       ch_faa_raw.map { _meta, faa -> faa },
       ch_map_file.map { _meta, map_file -> map_file },
       'plasmid'
    )

    ch_virus_prophage_fna = SEPARATE_VIRAL_SEQUENCES.out.map{ r -> r.fna }
        .mix( SEPARATE_PROPHAGES.out.map{ r -> r.fna } )
        .collectFile( name: 'viruses_and_prophages.fna' )
        .map{ fna -> [[id: 'combined'], fna] }

    //
    // ----------- Evaluate quality and rRNA content for viruses + prophages only
    // (plasmids are not checked for viral quality/rRNA contamination)
    //
    if ( ! params.skip_checkv ) {
        CHUNK_FNA (
            ch_virus_prophage_fna,
            [],                                        // length: (disabled) max number of nucleotides per chunk
            params.nucleotide_fasta_chunksize,         // size: max number of sequences per chunk
        )
        def ch_fna_chunks = CHUNK_FNA.out.chunked_output.transpose()

        CHECKV_ENDTOEND (
            ch_fna_chunks,
            params.checkv_db
        )

        CONCATENATE_CHECKV (
            CHECKV_ENDTOEND.out.quality_summary.groupTuple(),
            'tsv',
            'tsv'
        )
        quality_ch = CONCATENATE_CHECKV.out.csv.map { _meta, quality -> quality }
    } else {
        quality_ch = []
    }

    //
    // ----- Detect rRNA sequences ------
    //
    if ( ! params.skip_rrna_detection ) {
        BARRNAP(
          ch_virus_prophage_fna.map { meta, fasta -> [meta, fasta, "bac"] }
        )

        rna_gff = BARRNAP.out.gff.map { _meta, gff -> gff }
    } else {
        rna_gff = []
    }

    //
    // ----------- Deduplicate across categories, leave unique and save metadata
    // Compares hashsums across viruses/prophages/plasmids together (in case the same
    // sequence was classified differently by different predictions), prioritising
    // metagenome over genome and MGnify over third-party.
    // Also filters out non-determined-quality viruses/prophages, and drops the
    // corresponding GFF/FAA records for any sequence that gets deduplicated away.
    //
    CHOOSE_SEQUENCES (
        ch_map_file.map { meta, map_file -> record(meta: meta, map_file: map_file) },
        SEPARATE_VIRAL_SEQUENCES.out,
        SEPARATE_PROPHAGES.out,
        SEPARATE_PLASMIDS.out,
        quality_ch,
        rna_gff
    )
    chosen = CHOOSE_SEQUENCES.out

    COMPRESS_VIRAL_SEQ_FILES (
        chosen
        .map { r ->
            [r.meta, [r.virus_fna, r.virus_gff, r.virus_faa]]
        }
        .transpose()
    )

    COMPRESS_PLASMIDS_FILES (
        chosen
        .map { r ->
            [r.meta, [r.plasmid_fna, r.plasmid_gff, r.plasmid_faa]]
        }
        .transpose()
    )

    COMPRESS_PROPHAGES_FILES (
        chosen
        .map { r ->
            [r.meta, [r.prophage_fna, r.prophage_gff, r.prophage_faa]]
        }
        .transpose()
    )

    emit:
    mapfile             = ch_map_file             // [id:combined, metadata.tsv]

    separated_prophages = SEPARATE_PROPHAGES.out.map { r -> [r.meta, r.fna, r.gff, r.faa] }  // [[id:combined_prophage], fna, gff, faa]
    separated_viruses   = SEPARATE_VIRAL_SEQUENCES.out.map { r -> [r.meta, r.fna, r.gff, r.faa] }  // [[id:combined_virus], fna, gff, faa]
    separated_plasmids  = SEPARATE_PLASMIDS.out.map { r -> [r.meta, r.fna, r.gff, r.faa] }  // [[id:combined_plasmid], fna, gff, faa]

    input_metadata      = chosen.map { r -> [r.meta, r.metadata] }           // [id:combined, combined_metadata.tsv]
    metadata            = chosen.map { r -> [r.meta, r.filtered_metadata] }  // [id:combined, combined_filtered.tsv]
    excluded_qc         = chosen.map { r -> [r.meta, r.excluded_metadata] }  // [id:combined, excluded_metadata.tsv]

    viral_sequences     = chosen.map { r -> [r.meta, r.virus_fna, r.virus_gff, r.virus_faa] }     // [id:combined, combined_virus_filtered.fna, gff, faa]
    prophages           = chosen.map { r -> [r.meta, r.prophage_fna, r.prophage_gff, r.prophage_faa] }  // [id:combined, combined_prophage_filtered.fna, gff, faa]
    plasmids            = chosen.map { r -> [r.meta, r.plasmid_fna, r.plasmid_gff, r.plasmid_faa] }   // [id:combined, combined_plasmid_filtered.fna, gff, faa]

    combined_gff        = chosen.map { r -> r.filtered_gff }

    combined_faa        = chosen.map { r -> r.filtered_faa }

    stats_json          = chosen.map { r -> [r.meta, r.stats] }

}
