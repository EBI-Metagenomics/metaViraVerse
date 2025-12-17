include { METACEREBERUS                    } from '../../modules/local/metacerberus'
include { PHAMMSEQS                        } from '../../modules/local/phammseqs'

include { PROTEIN_PREDICTION               } from './protein_prediction'


workflow PROTEINS_PROCESSING {

    take:
    sequences

    main:

    ch_versions = Channel.empty()

    PROTEIN_PREDICTION(
       sequences
    )
    ch_versions = ch_versions.mix(PROTEIN_PREDICTION.out.versions)

    PHAMMSEQS(
       PROTEIN_PREDICTION.out.proteins
    )
    ch_versions = ch_versions.mix(PHAMMSEQS.out.versions)
    // TODO: maybe run PHAMCLUST on small dataset

    METACEREBERUS(
       sequences
    )
    ch_versions = ch_versions.mix(METACEREBERUS.out.versions)

    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}