/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FIND_CONCATENATE as CONCATENATE_VITAP        } from '../../modules/nf-core/find/concatenate'
include { GENOMAD_ENDTOEND                             } from '../../modules/nf-core/genomad/endtoend'
include { GENERATE_TAXONOMY_TABLE as TAX_VIPHOGS       } from '../../modules/local/generate_taxonomy_table'
include { GENERATE_TAXONOMY_TABLE as TAX_VITAP         } from '../../modules/local/generate_taxonomy_table'
include { VITAP                                        } from '../../modules/local/vitap'

include { VIPHOGS_ANNOTATION                           } from './viphogs_annotate'
include { TAXONOMY_VISUALISATION as VIS_VIPHOGS        } from './taxonomy_visualisation'
include { TAXONOMY_VISUALISATION as VIS_VITAP          } from './taxonomy_visualisation'



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow TAXONOMY_ASSIGNMENT {

    take:
    reps_stats_tsv
    combined_metadata
    ch_fna_chunks
    reps_faa
    reps_gff
    skip_vitap
    skip_genomad
    skip_viphogs
    viphog_db
    additional_model_data
    ncbi_db
    factor_file

    main:

    ch_versions = channel.empty()
    vitap_best  = channel.empty()

    //
    // Taxonomy genoMAD
    //
    if ( !skip_genomad ) {
        GENOMAD_ENDTOEND (
            ch_fna_chunks,
            params.genomad_db
        )
    }

    //
    // Taxonomy Viphogs
    //
    if ( !skip_viphogs ) {
        VIPHOGS_ANNOTATION (
            reps_faa,
            reps_gff,
            viphog_db,
            additional_model_data,
            ncbi_db,
            factor_file
        )
    }

    TAX_VIPHOGS (
        reps_stats_tsv
        .map { meta, table ->
            def new_meta = meta.clone()
            new_meta.tool = 'viphogs'
            tuple(new_meta, table)
        },
        combined_metadata
    )
    ch_versions = ch_versions.mix(TAX_VIPHOGS.out.versions)

    VIS_VIPHOGS(
       TAX_VIPHOGS.out.taxonomy_and_metadata
    )
    ch_versions = ch_versions.mix(VIS_VIPHOGS.out.versions)

    if ( !skip_vitap ) {
        //
        // Taxonomy VITAP
        //
        VITAP (
            ch_fna_chunks,
            params.vitap_db
        )
        ch_versions = ch_versions.mix(VITAP.out.versions)

        CONCATENATE_VITAP (
            VITAP.out.best_lineages.groupTuple(),
            1
        )

        TAX_VITAP (
            CONCATENATE_VITAP.out.file_out.map { meta, table ->
                def new_meta = meta.clone()
                new_meta.tool = 'vitap'
                tuple(new_meta, table)
            },
            combined_metadata
        )
        ch_versions = ch_versions.mix(TAX_VITAP.out.versions)

        VIS_VITAP (
           TAX_VITAP.out.taxonomy_and_metadata
        )
        ch_versions = ch_versions.mix(VIS_VITAP.out.versions)
    }


    emit:
    vitap_best     = CONCATENATE_VITAP.out.file_out
    versions       = ch_versions
}
