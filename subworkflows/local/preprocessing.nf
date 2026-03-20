include { CHOOSE_SEQUENCES                                } from '../../modules/local/choose_sequences'
include { RENAME_CONTIGS                                 } from '../../modules/local/rename_contigs'
include { SEPARATE_SEQUENCES as SEPARATE_VIRAL_SEQUENCES } from '../../modules/local/separate_sequences'
include { SEPARATE_SEQUENCES as SEPARATE_PLASMIDS        } from '../../modules/local/separate_sequences'
include { SEPARATE_SEQUENCES as SEPARATE_PROPHAGES       } from '../../modules/local/separate_sequences'

include { BARRNAP                                        } from '../../modules/nf-core/barrnap'

workflow PREPROCESSING {

    // TODO: process separately ASA, MAGs and third party and then combine

    take:
    input

    main:
    ch_versions = channel.empty()

    //
    // --- Aggregate catalogue MAGs and ASA results
    // Deduplicate sequences across samples, prioritising assembly over MAG
    //
    ch_fna_files = input.map { meta, gff, fna, faa, type, biome -> fna }.collect()
    ch_types     = input.map { meta, gff, fna, faa, type, biome -> type }.collect()
    ch_biomes    = input.map { meta, gff, fna, faa, type, biome -> biome }.collect()

    CHOOSE_SEQUENCES(
        ch_fna_files,
        ch_types,
        ch_biomes
    )
    ch_versions = ch_versions.mix(CHOOSE_SEQUENCES.out.versions)

    ch_fna_sequences = CHOOSE_SEQUENCES.out.combined_fna.map { fna -> tuple([id:'combined'], fna) }

    // TODO: implement that step for third party rename and probably ASA results
    //RENAME_CONTIGS(
    //   input
    //)

    //
    // ----- Detect rRNA sequences ------
    // TODO: implement that step for third party rename and probably ASA results
    if ( params.filter_rrna ) {
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

    all_viral_sequences = SEPARATE_VIRAL_SEQUENCES.out.chosen_sequences
        .map{ meta, seqs -> seqs }
        .combine(SEPARATE_PROPHAGES.out.chosen_sequences
        .map{ meta, seqs -> seqs })
        .flatMap { tuple -> tuple }
        .collectFile(name: "viral_sequences.fasta")
        .map{seqs -> [[id: 'viral_sequences'], seqs]}

    // publish
    all_viral_sequences.subscribe{ meta, seqs ->
            seqs.copyTo("${params.outdir}/${meta.id}/viral_sequences.fasta")
        }

    all_plasmids = SEPARATE_PLASMIDS.out.chosen_sequences
       .map{ meta, seqs -> seqs }
       .collectFile(name: "plasmids.fasta")
       .map{seqs -> [[id: 'plasmids'], seqs]}

    // publish
    all_plasmids.subscribe{ meta, seqs ->
            seqs.copyTo("${params.outdir}/${meta.id}/plasmids.fasta")
        }

    //all_gff = RENAME_CONTIGS.out.gff_renamed
    //   .map{ meta, gff -> gff }.collectFile(name: "combined.gff").map{gff -> [[id: 'combined'], gff]}

    //all_mapping = RENAME_CONTIGS.out.map_file
    //   .map{ meta, mapfile -> mapfile }.collectFile(name: "combined.map").map{mapfile -> [[id: 'combined'], mapfile]}

    emit:
    metadata       = CHOOSE_SEQUENCES.out.metadata
                     .map{ metadata -> [[id: 'combined'], metadata] }      // channel: [ [id: 'combined'], combined_meta.tsv ]
    viral_seqs     = all_viral_sequences
    plasmids       = all_plasmids
    all_gff        = input.map { meta, gff, fna, faa, type, biome -> gff }.collectFile( name: 'combined.gff' )
    all_faa        = input.map { meta, gff, fna, faa, type, biome -> faa }.collectFile( name: 'combined.faa' )
    all_mapping    = Channel.empty()
    versions       = ch_versions                        // channel: [ path(versions.yml) ]
}
