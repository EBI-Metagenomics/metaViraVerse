include { CHOOSE_SEQUENCES                               } from '../../modules/local/choose_sequences'
include { RENAME_CONTIGS as RENAME_CONTIGS_TMP           } from '../../modules/local/rename_contigs'
include { SEPARATE_SEQUENCES as SEPARATE_VIRAL_SEQUENCES } from '../../modules/local/separate_sequences'
include { SEPARATE_SEQUENCES as SEPARATE_PLASMIDS        } from '../../modules/local/separate_sequences'
include { SEPARATE_SEQUENCES as SEPARATE_PROPHAGES       } from '../../modules/local/separate_sequences'

include { BARRNAP                                        } from '../../modules/nf-core/barrnap'
include { CHECKV_ENDTOEND                                } from '../../modules/nf-core/checkv/endtoend'
include { THIRD_PARTY_DATA                               } from './third_party_data'

workflow PREPROCESSING {

    take:
    input
    third_party_input

    main:
    ch_versions = channel.empty()

    //
    // ----------- Preprocess third party data (coming not from MGnify)
    //
    THIRD_PARTY_DATA (
         third_party_input 
    )

    ch_fna_files = input.map { meta, gff, fna, faa -> tuple([meta, fna])}
    map_file = channel.empty()

    // TODO add third party data outputs
    // Review renaming for third party
    // review rename for ASA inputs
    //
    // ----------- Assign unique identifiers to all coming sequences
    //
    if ( !params.skip_rename ) {
        rename_input = input.map { meta, gff, fna, faa -> tuple([meta, fna, gff]) }
        RENAME_CONTIGS_TMP(
            rename_input,
            params.start_accession,
            params.end_accession
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
        map_file
    )

    ch_versions = ch_versions.mix(CHOOSE_SEQUENCES.out.versions)

    ch_fna_sequences = CHOOSE_SEQUENCES.out.filtered_fna

    //
    // ----- SEPARATE SEQUENCES INTO VIRAL AND PLASMIDS ------
    //
    SEPARATE_VIRAL_SEQUENCES(
       ch_fna_sequences,
       "viral_sequence"
    )
    ch_versions = ch_versions.mix(SEPARATE_VIRAL_SEQUENCES.out.versions)

    SEPARATE_PROPHAGES(
       ch_fna_sequences,
       "prophage"
    )

    SEPARATE_PLASMIDS(
       ch_fna_sequences,
       "plasmid"
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
