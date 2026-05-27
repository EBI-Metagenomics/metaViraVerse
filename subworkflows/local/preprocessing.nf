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

    ch_fna_files = input.map { meta, gff, fna, faa -> fna}

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
        ch_fna_files = RENAME_CONTIGS_TMP.out.contigs_renamed.map{meta, fna -> fna}
    }

    combined_fna = ch_fna_files
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
    ch_grouped = ch_fna_files
        .collect()               // produces one list: [[meta1, fna1], [meta2, fna2], ...]
        .map { tuples ->
            def fnas   = tuples.collect { it[1] }
            def types  = tuples.collect { it[0].type }
            def biomes = tuples.collect { it[0].biome }
            [fnas, types, biomes]
        }

    ch_types     = ch_grouped.map { fnas, types, biomes -> types }
    ch_biomes    = ch_grouped.map { fnas, types, biomes -> biomes }

    CHOOSE_SEQUENCES (
        ch_fna_files.collect(),
        CHECKV_ENDTOEND.out.quality_summary,
        ch_types,
        ch_biomes,

    )
    ch_versions = ch_versions.mix(CHOOSE_SEQUENCES.out.versions)

    ch_fna_sequences = CHOOSE_SEQUENCES.out.combined_fna.map { fna -> tuple([id:'combined'], fna) }

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

    combined_gff     = input.map { _meta, gff, _fna, _faa -> gff }
                        .collectFile( name: 'combined.gff' )
                        //.mix( THIRD_PARTY_DATA.out.gff.map { meta, gff -> gff } )

    combined_faa     = input.map { _meta, _gff, _fna, faa -> faa }
                        .collectFile( name: 'combined.faa' )
                        //.mix( THIRD_PARTY_DATA.out.faa.map { meta, faa -> faa } )

    versions         = ch_versions                        // channel: [ path(versions.yml) ]
}
