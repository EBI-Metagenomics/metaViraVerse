include { RENAME_CONTIGS                                 } from '../../modules/local/rename_contigs'
include { SEPARATE_SEQUENCES as SEPARATE_VIRAL_SEQUENCES } from '../../modules/local/separate_sequences'
include { SEPARATE_SEQUENCES as SEPARATE_PLASMIDS        } from '../../modules/local/separate_sequences'

workflow PREPROCESSING {

    take:
    input

    main:
    ch_versions = Channel.empty()

    RENAME_CONTIGS(
       input
    )

    SEPARATE_VIRAL_SEQUENCES(
       RENAME_CONTIGS.out.contigs_renamed,
       "viral_sequence"
    )

    SEPARATE_PLASMIDS(
       RENAME_CONTIGS.out.contigs_renamed,
       "plasmid"
    )

    all_viral_sequences = SEPARATE_VIRAL_SEQUENCES.out.chosen_sequences
       .map{ meta, seqs -> seqs }.collectFile(name: "viral_sequences.fasta").map{seqs -> [[id: 'viral_sequences'], seqs]}

    all_plasmids = SEPARATE_PLASMIDS.out.chosen_sequences
       .map{ meta, seqs -> seqs }.collectFile(name: "plasmids.fasta").map{seqs -> [[id: 'plasmids'], seqs]}

    all_gff = RENAME_CONTIGS.out.gff_renamed
       .map{ meta, gff -> gff }.collectFile(name: "combined.gff").map{gff -> [[id: 'combined'], gff]}

    all_mapping = RENAME_CONTIGS.out.map_file
       .map{ meta, mapfile -> mapfile }.collectFile(name: "combined.map").map{mapfile -> [[id: 'combined'], mapfile]}

    emit:
    viral_seqs     = all_viral_sequences
    plasmids       = all_plasmids
    all_gff        = all_gff
    all_mapping    = all_mapping
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}