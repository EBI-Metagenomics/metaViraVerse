include { METACEREBERUS                    } from '../../modules/local/metacerberus'
include { PHAMMSEQS                        } from '../../modules/local/phammseqs'

include { AMR_ANNOTATION                   } from '../ebi-metagenomics/amr_annotation'


workflow PROTEINS_PROCESSING {

    take:
    input  // (meta, reps_faa.uncompressed, reps_gff)

    main:

    ch_versions = channel.empty()

    ch_proteins = input.map{ meta, faa, _gff -> tuple(meta, faa) }

    //
    // -------- Antimicrobial resistence detection
    //
    AMR_ANNOTATION (
        input,
        params.amrfinderplus_db,
        params.deeparg_db,
        params.deeparg_db_version,
        params.deeparg_model,
        params.deeparg_tool_version,
        params.rgi_db,
        params.skip_amrfinderplus,
        params.skip_deeparg,
        params.skip_rgi
    )

    if (params.phammseqs) {
        //
        // -------- Assort phage protein sequences into phamilies using MMseqs2
        //
        PHAMMSEQS(
           ch_proteins
        )
        ch_versions = ch_versions.mix(PHAMMSEQS.out.versions)
        // TODO: maybe run PHAMCLUST on small dataset
    }

    //
    // -------- Functional annotation
    //
    METACEREBERUS(
       ch_proteins,
       params.metacerberus_db
    )
    ch_versions = ch_versions.mix(METACEREBERUS.out.versions)

    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}