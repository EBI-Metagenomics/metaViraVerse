include { CHOOSE_SEQUENCES                               } from '../../modules/local/choose_sequences'
include { RENAME_CONTIGS as RENAME_CONTIGS_TMP           } from '../../modules/local/rename_contigs'
include { RENAME_CONTIGS as RENAME_CONTIGS_COMBINED      } from '../../modules/local/rename_contigs'
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

    ch_fna_files = input.map { meta, gff, fna, faa -> tuple([meta, fna])}
    map_file = channel.empty()

    //
    // Review renaming for third party
    // review rename for ASA inputs
    //
    // ----------- Assign unique identifiers to all coming sequences by each record in samplesheet
    // We can't combine all sequences together here because we need to record it's biome and type later
    // That is why we have to give sequences temporary names (checkV can complain on original)
    // Those temporary names would be {meta.id}, ex. >barley1
    // No need to rename GFF here because GFF file would not be used in pre-processing
    //
    if ( !params.skip_rename ) {
        rename_input = input.map { meta, gff, fna, faa -> tuple([meta, fna, null]) }
        RENAME_CONTIGS_TMP(
            rename_input,
            false,
            false
        )
        ch_fna_files = RENAME_CONTIGS_TMP.out.contigs_renamed
        map_file = RENAME_CONTIGS_TMP.out.map_file
    }

    combined_fna = ch_fna_files.map{meta, fna -> fna}
        .collectFile(name: "input.fna")
        .map { fna -> tuple([id: 'combined'], fna) }

    //
    // ----------- Evaluate a quality for all coming sequences
    //
    CHECKV_ENDTOEND (
        combined_fna,
        params.checkv_db
    )

    //
    // ----- Detect rRNA sequences ------
    //
    if ( ! params.skip_rrna_detection ) {
        BARRNAP(
          combined_fna.map {id, fasta -> [id, fasta, "bac"]}
        )
        ch_versions = ch_versions.mix(BARRNAP.out.versions)

        rna_gff = BARRNAP.out.gff
    } else {
        rna_gff = tuple([id:'combined'], [])
    }

    //
    // ----------- Filter sequences, leave unique and save metadata
    // Deduplicate sequences across samples, prioritising assembly over MAG
    //
    ch_types     = ch_fna_files.map { meta, fna -> tuple([meta.type]) }
    ch_biomes    = ch_fna_files.map { meta, fna -> tuple([meta.biome]) }

    CHOOSE_SEQUENCES (
        ch_fna_files.map { meta, fna -> fna }.collect(),
        CHECKV_ENDTOEND.out.quality_summary,
        ch_types,
        ch_biomes,
        rna_gff,
        map_file.map { meta, fna -> fna }.collect()
    )

    ch_versions = ch_versions.mix(CHOOSE_SEQUENCES.out.versions)

    //
    // ----------- Rename all filtered sequences into one space
    // On that step we have all unique sequences passed quality control
    // They should now have unique identifiers coming from params.rename_accession
    //
    combined_gff = input.map { meta, gff, fna, faa -> gff }
        .collectFile(name: "input.gff")
        .map { gff -> tuple([id: 'combined'], gff) }
    RENAME_CONTIGS_COMBINED (
        CHOOSE_SEQUENCES.out.filtered_fna.join(combined_gff),
        params.start_accession,
        params.end_accession
    )

    // TODO add two mapping files combined file with original, temporary and unique
    ch_fna_sequences = RENAME_CONTIGS_COMBINED.out.contigs_renamed
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
