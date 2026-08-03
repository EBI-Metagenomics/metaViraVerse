include { CHOOSE_SEQUENCES                               } from '../../modules/local/choose_sequences'
include { RENAME_CONTIGS                                 } from '../../modules/local/rename_contigs'
include { SEPARATE_SEQUENCES as SEPARATE_VIRAL_SEQUENCES } from '../../modules/local/separate_sequences'
include { SEPARATE_SEQUENCES as SEPARATE_PLASMIDS        } from '../../modules/local/separate_sequences'
include { SEPARATE_SEQUENCES as SEPARATE_PROPHAGES       } from '../../modules/local/separate_sequences'

include { BARRNAP                                        } from '../../modules/nf-core/barrnap'
include { CSVTK_CONCAT as CONCATENATE_CHECKV             } from '../../modules/nf-core/csvtk/concat'
include { CHECKV_ENDTOEND                                } from '../../modules/nf-core/checkv/endtoend'
include { SEQKIT_SPLIT2 as CHUNK_FNA                     } from '../../modules/nf-core/seqkit/split2'


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

    ch_fna       = input.map { meta, gff, fna, faa -> fna }.collect()
    ch_gff       = input.map { meta, gff, fna, faa -> gff }.collect()
    ch_types     = input.map { meta, gff, fna, faa -> tuple([meta.source]) }.collect()
    ch_biomes    = input.map { meta, gff, fna, faa -> tuple([meta.biome]) }.collect()

    ch_fna_tp    = ch_third_party_data.map { meta, fna, gff, faa -> fna }.collect()
    ch_gff_tp    = ch_third_party_data.map { meta, fna, gff, faa -> gff }.collect()
    ch_types_tp  = ch_third_party_data.map { meta, fna, gff, faa -> tuple([meta.type]) }.collect()
    ch_biomes_tp = ch_third_party_data.map { meta, fna, gff, faa -> tuple([meta.biome]) }.collect()
    ch_source_tp = ch_third_party_data.map { meta, fna, gff, faa -> tuple([meta.source]) }.collect()

    // Raw pool of protein sequences (original, unrenamed protein IDs) from both sources,
    // used later to pull out just the proteins that survive separation + deduplication
    ch_faa_raw = input.map { meta, gff, fna, faa -> faa }
        .mix( ch_third_party_data.map { meta, fna, gff, faa -> faa } )
        .collectFile( name: 'combined_raw.faa' )
        .map { faa -> [[id: 'combined'], faa] }

    // TODO review how to handle meta
    RENAME_CONTIGS(
        ch_fna,
        ch_gff,
        ch_types,
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
    ch_versions = ch_versions.mix(RENAME_CONTIGS.out.versions)
    mapping = RENAME_CONTIGS.out.map_file.map{ map -> [[id: 'combined'], map] }

    ch_combined_fna = RENAME_CONTIGS.out.fna_renamed.map{ fna -> [[id: 'combined'], fna] }
    ch_combined_gff = RENAME_CONTIGS.out.gff_renamed.map{ gff -> [[id: 'combined'], gff] }

    //
    // ----- SEPARATE SEQUENCES INTO VIRUSES, PROPHAGES AND PLASMIDS ------
    // Splitting happens right after renaming, based on the 'definition' column of the
    // rename map (third_party_* definitions are folded into the same bucket as their
    // MGnify counterpart). FNA, GFF and FAA are all split together so downstream steps
    // (quality evaluation, choose_sequences) have a consistent triple per category.
    //

    SEPARATE_VIRAL_SEQUENCES(
       ch_combined_fna,
       ch_combined_gff,
       ch_faa_raw,
       mapping,
       'virus'
    )
    ch_versions = ch_versions.mix(SEPARATE_VIRAL_SEQUENCES.out.versions)

    SEPARATE_PROPHAGES(
       ch_combined_fna,
       ch_combined_gff,
       ch_faa_raw,
       mapping,
       'prophage'
    )
    ch_versions = ch_versions.mix(SEPARATE_PROPHAGES.out.versions)

    SEPARATE_PLASMIDS(
       ch_combined_fna,
       ch_combined_gff,
       ch_faa_raw,
       mapping,
       'plasmid'
    )
    ch_versions = ch_versions.mix(SEPARATE_PLASMIDS.out.versions)

    //
    // ----------- Evaluate quality and rRNA content for viruses + prophages only
    // (plasmids are not checked for viral quality/rRNA contamination)
    //

    ch_virus_prophage_fna = SEPARATE_VIRAL_SEQUENCES.out.chosen_sequences.map{ meta, fna -> fna }
        .mix( SEPARATE_PROPHAGES.out.chosen_sequences.map{ meta, fna -> fna } )
        .collectFile( name: 'viruses_and_prophages.fna' )
        .map{ fna -> [[id: 'combined'], fna] }

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
        SEPARATE_VIRAL_SEQUENCES.out.chosen_gff,
        SEPARATE_VIRAL_SEQUENCES.out.chosen_faa,
        SEPARATE_PROPHAGES.out.chosen_sequences,
        SEPARATE_PROPHAGES.out.chosen_gff,
        SEPARATE_PROPHAGES.out.chosen_faa,
        SEPARATE_PLASMIDS.out.chosen_sequences,
        SEPARATE_PLASMIDS.out.chosen_gff,
        SEPARATE_PLASMIDS.out.chosen_faa,
        CONCATENATE_CHECKV.out.csv,
        rna_gff,
        mapping
    )
    ch_versions = ch_versions.mix(CHOOSE_SEQUENCES.out.versions)

    emit:
    input_metadata   = CHOOSE_SEQUENCES.out.metadata           // [id:combined, combined_metadata.tsv]
    metadata         = CHOOSE_SEQUENCES.out.filtered_metadata  // [id:combined, combined_filtered.tsv]
    excluded_qc      = CHOOSE_SEQUENCES.out.excluded_metadata  // [id:combined, excluded_metadata.tsv]
    mapfile          = mapping                                 // [id:combined, metadata.tsv]

    viral_sequences  = CHOOSE_SEQUENCES.out.virus_fna     // [id:combined, combined_virus_filtered.fna]
    prophages        = CHOOSE_SEQUENCES.out.prophage_fna  // [id:combined, combined_prophage_filtered.fna]
    plasmids         = CHOOSE_SEQUENCES.out.plasmid_fna   // [id:combined, combined_plasmid_filtered.fna]

    combined_gff     = CHOOSE_SEQUENCES.out.filtered_gff
                        .map { meta, gff -> gff }

    combined_faa     = CHOOSE_SEQUENCES.out.filtered_faa
                        .map { meta, faa -> faa }

    versions         = ch_versions                        // channel: [ path(versions.yml) ]
}
