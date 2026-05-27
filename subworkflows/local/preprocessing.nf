include { CHOOSE_SEQUENCES                               } from '../../modules/local/choose_sequences'
include { RENAME_CONTIGS                                 } from '../../modules/local/rename_contigs'
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

    //
    // ----------- Assign unique identifiers to all coming sequences
    //
    if ( !params.skip_rename ) {
        rename_input = input.map { meta, gff, fna, faa, type, biome -> tuple([meta, fna, gff]) }
            .combine( THIRD_PARTY_DATA.out.fna.join(THIRD_PARTY_DATA.out.gff) )
        rename_input.view()
        RENAME_CONTIGS(
            rename_input,
            params.rename_accession,
            params.start_accession,
            params.end_accession
        )
    }

    //
    // ----------- Evaluate a quality for all coming sequences
    //
    ch_fna_files = input.map { meta, gff, fna, faa, type, biome -> fna}.mix(THIRD_PARTY_DATA.out.fna)
    CHECKV_ENDTOEND (
        ch_fna_files
            .collectFile(name: "input.fna")
            .map { fna -> tuple([id: 'combined'], fna) },
        params.checkv_db
    )

    //
    // ----------- Filter sequences, leave unique and save metadata
    // Deduplicate sequences across samples, prioritising assembly over MAG
    //
    ch_types     = input.map { meta, gff, fna, faa, type, biome -> type }
        .mix( THIRD_PARTY_DATA.out.types )
        .collect()
    ch_biomes    = input.map { meta, gff, fna, faa, type, biome -> biome }
        .mix( THIRD_PARTY_DATA.out.biomes )
        .collect()

    CHOOSE_SEQUENCES (
        ch_fna_files.collect(),
        ch_types,
        ch_biomes
    )
    ch_versions = ch_versions.mix(CHOOSE_SEQUENCES.out.versions)

    ch_fna_sequences = CHOOSE_SEQUENCES.out.combined_fna.map { fna -> tuple([id:'combined'], fna) }

    //
    // ----- Detect rRNA sequences ------
    // TODO: implement that step for third party rename and probably ASA results
    if ( ! params.skip_rrna_detection ) {
        BARRNAP(
          ch_fna_sequences.map {id, fasta -> [id, fasta, "bac"]}
        )
        ch_versions = ch_versions.mix(BARRNAP.out.versions)

        rna_gff = BARRNAP.out.gff
    } else {
        rna_gff = tuple([id:'combined'], [])
    }

    //
    // ----- SEPARATE SEQUENCES INTO VIRAL AND PLASMIDS ------
    //
    SEPARATE_VIRAL_SEQUENCES(
       CHOOSE_SEQUENCES.out.combined_fna.map { fna -> tuple([id:'combined'], fna) },
       "viral_sequence",
       rna_gff
    )
    ch_versions = ch_versions.mix(SEPARATE_VIRAL_SEQUENCES.out.versions)

    SEPARATE_PROPHAGES(
       CHOOSE_SEQUENCES.out.combined_fna.map { fna -> tuple([id:'combined'], fna) },
       "prophage",
       rna_gff
    )

    SEPARATE_PLASMIDS(
       CHOOSE_SEQUENCES.out.combined_fna.map { fna -> tuple([id:'combined'], fna) },
       "plasmid",
       rna_gff
    )
    ch_versions = ch_versions.mix(SEPARATE_PLASMIDS.out.versions)


    emit:
    metadata         = CHOOSE_SEQUENCES.out.metadata
                         .map{ metadata -> [[id: 'combined'], metadata] }      // channel: [ [id: 'combined'], combined_meta.tsv ]

    viral_sequences  = SEPARATE_VIRAL_SEQUENCES.out.chosen_sequences
    prophages        = SEPARATE_PROPHAGES.out.chosen_sequences
    plasmids         = SEPARATE_PLASMIDS.out.chosen_sequences

    combined_gff     = input.map { _meta, gff, _fna, _faa, _type, _biome -> gff }
                        .mix( THIRD_PARTY_DATA.out.gff.map { meta, gff -> gff } )
                        .collectFile( name: 'combined.gff' )
    combined_faa     = input.map { _meta, _gff, _fna, faa, _type, _biome -> faa }
                        .mix( THIRD_PARTY_DATA.out.faa.map { meta, faa -> faa } )
                        .collectFile( name: 'combined.faa' )

    versions         = ch_versions                        // channel: [ path(versions.yml) ]
}
