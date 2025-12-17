include { COMBINED_GENE_CALLER    } from '../../modules/local/cgc'
include { PHANOTATE               } from '../../modules/local/phanotate'

include { PRODIGAL                } from '../../modules/nf-core/prodigal'

workflow PROTEIN_PREDICTION {

    take:
    sequences

    main:

    ch_versions = Channel.empty()

    PRODIGAL(
       sequences,
       "gff"
    )
    ch_versions = ch_versions.mix(PRODIGAL.out.versions)

    PHANOTATE(
       sequences
    )
    ch_versions = ch_versions.mix(PHANOTATE.out.versions)

    COMBINED_GENE_CALLER(
       PRODIGAL.out.amino_acid_fasta,
       PHANOTATE.out.proteins
    )
    ch_versions = ch_versions.mix(COMBINED_GENE_CALLER.out.versions)

    emit:
    proteins       = COMBINED_GENE_CALLER.out.combined_proteins
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}