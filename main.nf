#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    EBI-Metagenomics/metaviraverse
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/EBI-Metagenomics/metaviraverse
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

feature/OncoProteinStructure
include { VIRAVERSE  } from './workflows/viraverse'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_viraverse_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_viraverse_pipeline'

include { PROTEIN_STRUCTURE } from './subworkflows/local/protein_structure/main'
=======
include { METAVIRAVERSE           } from './workflows/metaviraverse'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_metaviraverse_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_metaviraverse_pipeline'
dev
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow EBIMETAGENOMICS_METAVIRAVERSE {

    take:
    samplesheet // channel: samplesheet read in from --input

    main:

    //
    // WORKFLOW: Run pipeline
    //
    METAVIRAVERSE (
        samplesheet
    )
    emit:
    multiqc_report = METAVIRAVERSE.out.multiqc_report // channel: /path/to/multiqc_report.html
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input
    )

    //
    // WORKFLOW: Run main workflow
    //
    EBIMETAGENOMICS_METAVIRAVERSE (
        PIPELINE_INITIALISATION.out.samplesheet
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.outdir,
        params.monochrome_logs,
        EBIMETAGENOMICS_METAVIRAVERSE.out.multiqc_report
    )
    //
    // OPTIONAL: Protein structure prediction and annotation
    // Activated when params.run_protein_structure = true
    // and .faa is provided in the samplesheet
    //
    if (params.run_protein_structure) {
        ch_faa = ch_samplesheet
            .map { meta, reads, faa -> faa ? [ meta, file(faa) ] : null }
            .filter { it != null }

        PROTEIN_STRUCTURE (
            ch_faa,
            ch_reps_stats,          // from upstream clustering process
            params.bfvd_db      ?: 'BFVD',
            params.plddt_cutoff ?: 0.7
        )
    }

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
