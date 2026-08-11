include { CHOOSE_SEQUENCES                               } from '../../../modules/local/choose_sequences'
include { RENAME_CONTIGS                                 } from '../../../modules/local/rename_contigs'
include { SEPARATE_SEQUENCES as SEPARATE_VIRAL_SEQUENCES } from '../../../modules/local/separate_sequences'
include { SEPARATE_SEQUENCES as SEPARATE_PLASMIDS        } from '../../../modules/local/separate_sequences'
include { SEPARATE_SEQUENCES as SEPARATE_PROPHAGES       } from '../../../modules/local/separate_sequences'

include { BARRNAP                                        } from '../../../modules/nf-core/barrnap'
include { CSVTK_CONCAT as CONCATENATE_CHECKV             } from '../../../modules/nf-core/csvtk/concat'
include { CHECKV_ENDTOEND                                } from '../../../modules/nf-core/checkv/endtoend'
include { SEQKIT_SPLIT2 as CHUNK_FNA                     } from '../../../modules/nf-core/seqkit/split2'


workflow PREPROCESSING {

    take:
    input   // [meta, gff, fna, faa]
    ch_third_party_data  // [meta[id, biome, type, source], fasta, gff, faa]

    main:
    ch_versions = channel.empty()

    // ----------- Assign unique identifiers to all coming sequences by each record in samplesheet
    // We will combine all sequences together and keep track of names, biomes and types in map file
    // The new names would be {prefix}{number}, ex. >seq1
    // For MGnify processing rename_accession should be MGYV
    // It will also rename ID in attributes column in GFF
    // That step is running with --combine option and it will return one renamed FASTA and GFF

    ch_fna       = input.map { meta, fna, gff, faa -> fna }.collect()
    ch_gff       = input.map { meta, fna, gff, faa -> gff }.collect()
    ch_sources   = input.map { meta, fna, gff, faa -> tuple([meta.source]) }.collect()
    ch_biomes    = input.map { meta, fna, gff, faa -> tuple([meta.biome]) }.collect()

    ch_fna_tp    = ch_third_party_data.map { meta, fna, gff, faa -> fna }.collect()
    ch_gff_tp    = ch_third_party_data.map { meta, fna, gff, faa -> gff }.collect()
    ch_types_tp  = ch_third_party_data.map { meta, fna, gff, faa -> tuple([meta.type]) }.collect()
    ch_biomes_tp = ch_third_party_data.map { meta, fna, gff, faa -> tuple([meta.biome]) }.collect()
    ch_source_tp = ch_third_party_data.map { meta, fna, gff, faa -> tuple([meta.source]) }.collect()

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
        params.plasmid_identifier
    )

    //
    // ----- SEPARATE SEQUENCES INTO VIRUSES, PROPHAGES AND PLASMIDS ------
    // Splitting happens right after renaming, based on the 'definition' column of the
    // rename map (third_party_* definitions are folded into the same bucket as their
    // MGnify counterpart). FNA, GFF and FAA are all split together so downstream steps
    // (quality evaluation, choose_sequences) have a consistent triple per category.
    //

    SEPARATE_VIRAL_SEQUENCES(
       RENAME_CONTIGS.out.fna_renamed,
       RENAME_CONTIGS.out.gff_renamed,
       ch_faa_raw,
       RENAME_CONTIGS.out.map_file,
       'virus'
    )

    SEPARATE_PROPHAGES(
       RENAME_CONTIGS.out.fna_renamed,
       RENAME_CONTIGS.out.gff_renamed,
       ch_faa_raw,
       RENAME_CONTIGS.out.map_file,
       'prophage'
    )

    SEPARATE_PLASMIDS(
       RENAME_CONTIGS.out.fna_renamed,
       RENAME_CONTIGS.out.gff_renamed,
       ch_faa_raw,
       RENAME_CONTIGS.out.map_file,
       'plasmid'
    )

    ch_virus_prophage_fna = SEPARATE_VIRAL_SEQUENCES.out.chosen_sequences.map{ meta, fna, gff, faa -> fna }
        .mix( SEPARATE_PROPHAGES.out.chosen_sequences.map{ meta, fna, gff, faa -> fna } )
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
        ch_versions = ch_versions.mix(CHUNK_FNA.out.versions)
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
        quality_ch = CONCATENATE_CHECKV.out.csv
    } else {
        quality_ch = tuple([id:'combined'], [])
    }

    //
    // ----- Detect rRNA sequences ------
    //
    if ( ! params.skip_rrna_detection ) {
        BARRNAP(
          ch_virus_prophage_fna.map { meta, fasta -> [meta, fasta, "bac"] }
        )
        ch_versions = ch_versions.mix(BARRNAP.out.versions)

        rna_gff = BARRNAP.out.gff
    } else {
        rna_gff = tuple([id:'combined'], [])
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
        SEPARATE_VIRAL_SEQUENCES.out.chosen_sequences,
        SEPARATE_PROPHAGES.out.chosen_sequences,
        SEPARATE_PLASMIDS.out.chosen_sequences,
        quality_ch,
        rna_gff,
        RENAME_CONTIGS.out.map_file
    )

    emit:
    mapfile             = RENAME_CONTIGS.out.map_file             // [id:combined, metadata.tsv]

    separated_prophages = SEPARATE_PROPHAGES.out.chosen_sequences  // [[id:combined_prophage], fna, gff, faa]
    separated_viruses   = SEPARATE_VIRAL_SEQUENCES.out.chosen_sequences  // [[id:combined_virus], fna, gff, faa]
    separated_plasmids  = SEPARATE_PLASMIDS.out.chosen_sequences  // [[id:combined_plasmid], fna, gff, faa]

    input_metadata      = CHOOSE_SEQUENCES.out.metadata           // [id:combined, combined_metadata.tsv]
    metadata            = CHOOSE_SEQUENCES.out.filtered_metadata  // [id:combined, combined_filtered.tsv]
    excluded_qc         = CHOOSE_SEQUENCES.out.excluded_metadata  // [id:combined, excluded_metadata.tsv]

    viral_sequences     = CHOOSE_SEQUENCES.out.virus_data     // [id:combined, combined_virus_filtered.fna, gff, faa]
    prophages           = CHOOSE_SEQUENCES.out.prophage_data  // [id:combined, combined_prophage_filtered.fna, gff, faa]
    plasmids            = CHOOSE_SEQUENCES.out.plasmid_data   // [id:combined, combined_plasmid_filtered.fna, gff, faa]

    combined_gff        = CHOOSE_SEQUENCES.out.filtered_data.map { meta, fna, gff, faa -> gff }

    combined_faa        = CHOOSE_SEQUENCES.out.filtered_data.map { meta, fna, gff, faa -> faa }

    versions            = ch_versions                        // channel: [ path(versions.yml) ]
}
