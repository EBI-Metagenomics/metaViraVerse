#!/usr/bin/env nextflow
nextflow.enable.types = true

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

include { METAVIRAVERSE ; MetaviraverseParams } from './workflows/metaviraverse'
include { PreprocessingParams     } from './subworkflows/local/preprocessing/main'
include { ViralParams             } from './subworkflows/local/process_viral_sequences'
include { PlasmidParams           } from './subworkflows/local/process_plasmids'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_metaviraverse_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_metaviraverse_pipeline'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PIPELINE PARAMETERS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Params used by the pipeline code. Params that only affect configuration
    (publishing, institutional profiles, reports) are declared in nextflow.config,
    which also keeps `outdir`, `rename_accession` and `monochrome_logs` for the config.
    Descriptions: nextflow_schema.json.
*/

params {
    // Input/output
    input: Path?
    third_party_input: Path?
    catalogues_metadata: Path
    outdir: String

    // Pre-processing
    viral_sequence_identifier: String = 'viral_sequence'
    prophage_identifier: String = 'prophage'
    plasmid_identifier: String = 'plasmid'
    skip_rrna_detection: Boolean = false
    rename_accession: String = 'seq'
    start_accession: Integer = 1
    end_accession: Integer = 1

    // Quality evaluation
    checkv_db: Path?
    skip_checkv: Boolean = false

    // Chunking files for parallelization (passed to --by-size of seqkit split2)
    nucleotide_fasta_chunksize: Integer = 25000
    protein_annotation_fasta_chunksize: Integer = 300000
    nucleotide_fasta_chunksize_iphop: Integer = 5000

    // Taxonomy
    vitap_db: Path?
    skip_vitap: Boolean = false
    genomad_db: Path?
    skip_genomad: Boolean = false
    viphog_db: Path?
    additional_viphog_model_data: Path?
    ncbi_db: Path?
    viphogs_factor_file: Path?
    skip_viphogs: Boolean = false
    viphogs_evalue: Float = 0.01

    // Clustering
    cluster_vclust: Boolean = true

    // AMR detection
    amrfinderplus_db: String?
    deeparg_db: String?
    deeparg_db_version: String?
    deeparg_model: String?
    deeparg_tool_version: String?
    rgi_db: String?
    skip_amrfinderplus_for_viruses: Boolean = true
    skip_deeparg_for_viruses: Boolean = false
    skip_rgi_for_viruses: Boolean = false
    skip_amrfinderplus_for_plasmids: Boolean = true
    skip_deeparg_for_plasmids: Boolean = false
    skip_rgi_for_plasmids: Boolean = false

    // HMM search
    annotation_db: Path?

    // Host detection
    iphop_db: Path?
    skip_iphop: Boolean = false
    predict_host_from_custom_spacers: Boolean = false
    custom_spacers_fasta: Path?
    custom_spacers_metadata: Path?
    custom_genomes_metadata: Path?

    // Protein clustering
    skip_phammseqs: Boolean = false
    skip_mmseqs: Boolean = false
    mmseqs_id_threshold: Float = 0.9
    mmseqs_cov_threshold: Float = 0.8

    // Plasmids
    repsearch_db: Path?
    repfilter_db: Path?
    rna_inc_db: Path?
    incsearch_db: Path?
    mobsearch_db: Path?
    mobsuite_db: Path?

    // MultiQC
    multiqc_config: Path?
    multiqc_logo: String?

    // Boilerplate
    monochrome_logs: Boolean = false
    version: Boolean = false
    validate_params: Boolean = true
}

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
    samplesheet             // channel: samplesheet read in from --input
    samplesheet_third_party // channel: samplesheet read in from --third_party_input
    opts: MetaviraverseParams
    preprocessing_opts: PreprocessingParams
    viral_opts: ViralParams
    plasmid_opts: PlasmidParams

    main:

    //
    // WORKFLOW: Run pipeline
    //
    multiqc_report = METAVIRAVERSE (
        samplesheet.map { meta, fna, gff, faa -> record(meta: meta, fna: fna, gff: gff, faa: faa) },
        samplesheet_third_party.map { meta, fasta -> record(meta: meta, fasta: fasta) },
        opts,
        preprocessing_opts,
        viral_opts,
        plasmid_opts
    )

    emit:
    multiqc_report = multiqc_report // channel: /path/to/multiqc_report.html
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
    init = PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.third_party_input
    )

    //
    // WORKFLOW: Run main workflow
    //
    multiqc_report = EBIMETAGENOMICS_METAVIRAVERSE (
        init.samplesheet,
        init.samplesheet_third_party,
        params,
        params,
        params,
        params
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.outdir,
        params.monochrome_logs,
        multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
