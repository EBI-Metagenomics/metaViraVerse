/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FIND_CONCATENATE as CONCATENATE_VITAP        } from '../../../modules/nf-core/find/concatenate'
include { GENOMAD_ENDTOEND                             } from '../../../modules/nf-core/genomad/endtoend'
include { CSVTK_CONCAT as CONCATENATE_GENOMAD          } from '../../../modules/nf-core/csvtk/concat'

include { GENERATE_TAXONOMY_TABLE as TAX_VIPHOGS       } from '../../../modules/local/generate_taxonomy_table'
include { GENERATE_TAXONOMY_TABLE as TAX_VITAP         } from '../../../modules/local/generate_taxonomy_table'
include { GENERATE_TAXONOMY_TABLE as TAX_GENOMAD       } from '../../../modules/local/generate_taxonomy_table'
include { VITAP                                        } from '../../../modules/local/vitap'

include { VIPHOGS_ANNOTATION                           } from '../viphogs_annotate/main'
include { TAXONOMY_VISUALISATION as VIS_VIPHOGS        } from '../taxonomy_visualisation/main'
include { TAXONOMY_VISUALISATION as VIS_VITAP          } from '../taxonomy_visualisation/main'
include { TAXONOMY_VISUALISATION as VIS_GENOMAD        } from '../taxonomy_visualisation/main'



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow TAXONOMY_ASSIGNMENT {

    take:
    combined_metadata
    ch_fna_chunks
    reps_faa
    reps_gff
    skip_vitap
    skip_genomad
    skip_viphogs
    viphog_db
    additional_viphog_model_data
    ncbi_db
    viphogs_factor_file

    main:

    ch_versions      = channel.empty()
    vitap_best       = channel.empty()
    genomad_taxonomy = channel.empty()
    viphogs_taxonomy = channel.empty()

    //
    // Taxonomy genoMAD
    //
    if ( !skip_genomad ) {
        GENOMAD_ENDTOEND (
            ch_fna_chunks,
            params.genomad_db
        )

        CONCATENATE_GENOMAD (
            GENOMAD_ENDTOEND.out.virus_summary.groupTuple(),
            'tsv',
            'tsv'
        )

        TAX_GENOMAD (
            CONCATENATE_GENOMAD.out.csv.map { meta, table ->
                def new_meta = meta.clone()
                new_meta.tool = 'genomad'
                tuple(new_meta, table)
            },
            combined_metadata
        )
        ch_versions = ch_versions.mix(TAX_GENOMAD.out.versions)

        VIS_GENOMAD (
            TAX_GENOMAD.out.taxonomy_and_metadata
        )
        ch_versions = ch_versions.mix(VIS_GENOMAD.out.versions)

        genomad_taxonomy = CONCATENATE_GENOMAD.out.csv
    }

    //
    // Taxonomy Viphogs
    //
    if ( !skip_viphogs ) {
        VIPHOGS_ANNOTATION (
            reps_faa,
            reps_gff,
            viphog_db,
            additional_viphog_model_data,
            ncbi_db,
            viphogs_factor_file
        )

        TAX_VIPHOGS (
            VIPHOGS_ANNOTATION.out.assignment
            .map { meta, table ->
                def new_meta = meta.clone()
                new_meta.tool = 'viphogs'
                tuple(new_meta, table)
            },
            combined_metadata
        )
        ch_versions = ch_versions.mix(TAX_VIPHOGS.out.versions)

        VIS_VIPHOGS (
            TAX_VIPHOGS.out.taxonomy_and_metadata
        )
        ch_versions = ch_versions.mix(VIS_VIPHOGS.out.versions)

        viphogs_taxonomy = VIPHOGS_ANNOTATION.out.assignment
    }

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

        vitap_best = CONCATENATE_VITAP.out.file_out
    }


    emit:
    vitap_best       = vitap_best
    genomad_taxonomy = genomad_taxonomy
    viphogs_taxonomy = viphogs_taxonomy
    versions         = ch_versions
}
