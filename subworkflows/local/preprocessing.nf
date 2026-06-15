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

    ch_fna       = input.map { meta, gff, fna, faa -> fna }.collect()
    ch_gff       = input.map { meta, gff, fna, faa -> gff }.collect()
    ch_types     = input.map { meta, gff, fna, faa -> tuple([meta.type]) }.collect()
    ch_biomes    = input.map { meta, gff, fna, faa -> tuple([meta.biome]) }.collect()

    // TODO review how to handle meta
    RENAME_CONTIGS(
        ch_fna,
        ch_gff,
        params.start_accession,
        params.end_accession,
        ch_types,
        ch_biomes
    )
    ch_versions = ch_versions.mix(RENAME_CONTIGS.out.versions)
    mapping = RENAME_CONTIGS.out.map_file.map{ map -> [[id: 'combined'], map] }

    //
    // ----------- Evaluate a quality for all coming sequences
    //

    CHUNK_FNA (
        RENAME_CONTIGS.out.fna_renamed.map{ fna -> [[id: 'combined'], fna] },
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
          RENAME_CONTIGS.out.fna_renamed.map {fasta -> [[id: 'combined'], fasta, "bac"]}
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
        RENAME_CONTIGS.out.fna_renamed.map{ fna -> [[id: 'combined'], fna] },
        RENAME_CONTIGS.out.gff_renamed.map{ gff -> [[id: 'combined'], gff] },
        CONCATENATE_CHECKV.out.file_out,
        rna_gff,
        mapping
    )
    ch_versions = ch_versions.mix(CHOOSE_SEQUENCES.out.versions)

    ch_fna_sequences = CHOOSE_SEQUENCES.out.filtered_fna

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
    metadata         = CHOOSE_SEQUENCES.out.filtered_metadata  // [id:combined, metadata.tsv]
    excluded_qc      = CHOOSE_SEQUENCES.out.excluded_metadata  // [id:combined, excluded_metadata.tsv]
    mapfile          = mapping                        // [id:combined, metadata.tsv]

    viral_sequences  = SEPARATE_VIRAL_SEQUENCES.out.chosen_sequences  // [id:combined, viruses.fasta]
    prophages        = SEPARATE_PROPHAGES.out.chosen_sequences        // [id:combined, prophages.fasta]
    plasmids         = SEPARATE_PLASMIDS.out.chosen_sequences         // [id:combined, plasmids.fasta]

    combined_gff     = CHOOSE_SEQUENCES.out.filtered_gff
                        .map { meta, gff -> gff }
                        //.mix( THIRD_PARTY_DATA.out.gff.map { meta, gff -> gff } )

    combined_faa     = input.map { _meta, _gff, _fna, faa -> faa }
                        .collectFile( name: 'combined.faa' )
                        //.mix( THIRD_PARTY_DATA.out.faa.map { meta, faa -> faa } )

    versions         = ch_versions                        // channel: [ path(versions.yml) ]
}
