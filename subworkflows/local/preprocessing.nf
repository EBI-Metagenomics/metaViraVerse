include { CHOOSE_SEQUENCES                               } from '../../modules/local/choose_sequences'
include { RENAME_CONTIGS                                 } from '../../modules/local/rename_contigs'
include { SEPARATE_SEQUENCES as SEPARATE_VIRAL_SEQUENCES } from '../../modules/local/separate_sequences'
include { SEPARATE_SEQUENCES as SEPARATE_PLASMIDS        } from '../../modules/local/separate_sequences'
include { SEPARATE_SEQUENCES as SEPARATE_PROPHAGES       } from '../../modules/local/separate_sequences'

include { BARRNAP                                        } from '../../modules/nf-core/barrnap'
include { CHECKV_ENDTOEND                                } from '../../modules/nf-core/checkv/endtoend'

workflow PREPROCESSING {

    take:
    input   // [meta, gff, fna, faa]

    main:
    ch_versions = channel.empty()

    //
    // TODO Review renaming for third party
    // TODO review rename for ASA inputs
    //
    // ----------- Assign unique identifiers to all coming sequences by each record in samplesheet
    // We will combine all sequences together and keep track of names, biomes and types in map file
    // The new names would be {prefix}{number}, ex. >seq1
    // For MGnify processing rename_accession should be MGYV
    // It will also rename ID in attributes column in GFF
    // That step is running with --combine option and it will return one renamed FASTA and GFF

    rename_input = input.map { meta, gff, fna, faa -> tuple([meta, fna, gff]) }
    ch_types     = ch_fna_files.map { meta, fna -> tuple([meta.type]) }
    ch_biomes    = ch_fna_files.map { meta, fna -> tuple([meta.biome]) }

    RENAME_CONTIGS(
        rename_input,
        params.start_accession,
        params.end_accession,
        ch_types,
        ch_biomes
    )

    //
    // ----------- Evaluate a quality for all coming sequences
    //
    CHECKV_ENDTOEND (
        RENAME_CONTIGS.out.fna_renamed,
        params.checkv_db
    )

    //
    // ----- Detect rRNA sequences ------
    //
    if ( ! params.skip_rrna_detection ) {
        BARRNAP(
          RENAME_CONTIGS.out.fna_renamed.map {id, fasta -> [id, fasta, "bac"]}
        )
        ch_versions = ch_versions.mix(BARRNAP.out.versions)

        rna_gff = BARRNAP.out.gff
    } else {
        rna_gff = tuple([id:'combined'], [])
    }

    //
    // ----------- Filter sequences, leave unique and save metadata
    // Deduplicate sequences across samples, prioritising assembly over MAG
    // Remove non-determined quality viruses
    //
    CHOOSE_SEQUENCES (
        RENAME_CONTIGS.out.fna_renamed,
        RENAME_CONTIGS.out.gff_renamed,
        CHECKV_ENDTOEND.out.quality_summary,
        rna_gff,
        RENAME_CONTIGS.out.map_file
    )
    ch_versions = ch_versions.mix(CHOOSE_SEQUENCES.out.versions)

    ch_fna_sequences = CHOOSE_SEQUENCES.out.filtered_fna
    mapping = RENAME_CONTIGS_COMBINED.out.map_file

    //
    // ----- SEPARATE SEQUENCES INTO VIRAL AND PLASMIDS ------
    //
    SEPARATE_VIRAL_SEQUENCES(
       ch_fna_sequences,
       "viral_sequence",
       mapping
    )
    ch_versions = ch_versions.mix(SEPARATE_VIRAL_SEQUENCES.out.versions)

    SEPARATE_PROPHAGES(
       ch_fna_sequences,
       "prophage",
       mapping
    )

    SEPARATE_PLASMIDS(
       ch_fna_sequences,
       "plasmid",
       mapping
    )
    ch_versions = ch_versions.mix(SEPARATE_PLASMIDS.out.versions)


    emit:
    metadata         = CHOOSE_SEQUENCES.out.metadata
                         .map{ metadata -> [[id: 'combined'], metadata] }      // channel: [ [id: 'combined'], combined_meta.tsv ]

    viral_sequences  = SEPARATE_VIRAL_SEQUENCES.out.chosen_sequences
    prophages        = SEPARATE_PROPHAGES.out.chosen_sequences
    plasmids         = SEPARATE_PLASMIDS.out.chosen_sequences

    combined_gff     = input.map { _meta, gff, _fna, _faa -> gff }
                        .collectFile( name: 'combined.gff' )
                        //.mix( THIRD_PARTY_DATA.out.gff.map { meta, gff -> gff } )

    combined_faa     = input.map { _meta, _gff, _fna, faa -> faa }
                        .collectFile( name: 'combined.faa' )
                        //.mix( THIRD_PARTY_DATA.out.faa.map { meta, faa -> faa } )

    versions         = ch_versions                        // channel: [ path(versions.yml) ]
}
