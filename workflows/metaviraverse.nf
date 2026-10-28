nextflow.enable.types = true

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap                      } from 'plugin/nf-schema'
include { paramsSummaryMultiqc                  } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { processVersionsFromYAML               } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { workflowVersionToYAML                 } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                } from '../subworkflows/local/utils_nfcore_metaviraverse_pipeline'

include { PREPROCESSING                         } from '../subworkflows/local/preprocessing/main'
include { PROCESS_VIRAL_SEQUENCES               } from '../subworkflows/local/process_viral_sequences'
include { PROCESS_PLASMIDS                      } from '../subworkflows/local/process_plasmids'
include { THIRD_PARTY_DATA                      } from '../subworkflows/local/third_party_data'


include { COLLECT_CATALOGUE_STATS               } from '../modules/local/collect_catalogue_stats'
include { COLLECT_METADATA                      } from '../modules/local/collect_metadata'
include { SANKEY_PLOT as SANKEY_PLASMIDS_HOST   } from '../modules/local/sankey_plot/main'
include { SANKEY_PLOT as SANKEY_VIRUSES_HOST    } from '../modules/local/sankey_plot/main'

include { MULTIQC                               } from '../modules/nf-core/multiqc'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

record Sample {
    meta: Map     // [id, source, biome]
    fna: Path
    gff: Path
    faa: Path
}

record ThirdPartyInput {
    meta: Map     // [id, type, source, biome, study_accession, sample_accession]
    fasta: Path
}

record ThirdPartySample {
    meta: Map
    fna: Path
    gff: Path
    faa: Path
}

workflow METAVIRAVERSE {

    take:
    samples: Channel<Sample>                       // samplesheet read in from --input
    third_party_input: Channel<ThirdPartyInput>    // samplesheet read in from --third_party_input

    main:

    //
    // ----------- Preprocess third party data (coming not from MGnify)
    //
    def third_party: Channel<ThirdPartySample> = channel.empty()
    if (params.third_party_input) {
        third_party = THIRD_PARTY_DATA(
            third_party_input,
            params.outdir
        )
    }

    //
    // Separate viral sequences and plasmids
    //
    preprocessed = PREPROCESSING (
       samples,
       third_party,
       record(
           start_accession: params.start_accession,
           end_accession: params.end_accession,
           viral_sequence_identifier: params.viral_sequence_identifier,
           prophage_identifier: params.prophage_identifier,
           plasmid_identifier: params.plasmid_identifier,
           rename_accession: params.rename_accession
       )
    )

    viruses = preprocessed
        .flatMap { r -> [r.virus_fna, r.prophage_fna] }
        .collectFile(name: "viruses.fasta")
        .map { seqs -> record(meta: [id: 'viruses'], fna: seqs) }

    //
    // Process viruses
    //
    viral_catalogue = PROCESS_VIRAL_SEQUENCES (
       viruses,
       preprocessed.map { r -> r.filtered_gff },   // concatenated GFF built from all inputs
       preprocessed.map { r -> r.filtered_faa },   // concatenated FAA built from all inputs
       preprocessed.map { r -> r.filtered_metadata },
       preprocessed.map { r -> r.map_file }
    )

    // Process plasmids filtered from input, those already compressed and published in PREPROCESSING

    plasmids = preprocessed
       .flatMap { r -> [r.plasmid_fna] }
       .collectFile(name: "plasmids.fasta")
       .map { seqs -> record(meta: [id: 'plasmids'], fna: seqs) }

    //
    // Process plasmids
    //
    plasmid_catalogue = PROCESS_PLASMIDS (
       plasmids,
       preprocessed.map { r -> r.filtered_faa },
       preprocessed.map { r -> r.filtered_gff }   // concatenated GFF built from all inputs
    )

    // one record with both catalogues and the preprocessing results
    catalogues = viral_catalogue
        .combine(plasmid_catalogue)
        .map { viral, plasmid -> record(viral: viral, plasmid: plasmid) }
        .combine(preprocessed: preprocessed)

    //
    // Collect stats for whole catalogue into JSON
    //
    COLLECT_CATALOGUE_STATS (
        catalogues.map { r ->
            record(
                json: r.preprocessed.stats,
                metadata: r.preprocessed.filtered_metadata,
                clusters_viruses: r.viral.reps_list,
                clusters_plasmids: r.plasmid.reps_list,
                excluded_metadata: r.preprocessed.excluded_metadata,
                initial_metadata: r.preprocessed.metadata,
                clusters_viruses_prot: r.viral.faa_gz,
                clusters_plasmids_prot: r.plasmid.faa_gz
            )
        }
    )

    //
    // Collect viral and plasmid metadata
    //
    final_metadata = COLLECT_METADATA (
        catalogues.map { r ->
            record(
                combined_meta: r.preprocessed.filtered_metadata,
                viral_meta: r.viral.meta,
                viral_clusters: r.viral.clusters_tsv,
                plasmid_meta: r.plasmid.meta,
                plasmid_clusters: r.plasmid.clusters_tsv,
                // optional annotations and host predictions: null when the step was skipped
                vitap_best: r.viral.vitap_best,
                viphogs_taxonomy: r.viral.viphogs,
                genomad: r.viral.genomad,
                map_file: r.preprocessed.map_file,
                iphop_genome: r.viral.iphop_genome,
                iphop_genus: r.viral.iphop_genus,
                spacepharer: r.viral.spacepharer
            )
        },
        file(params.catalogues_metadata)
    )


    // Plot Sankey for hosts
    SANKEY_PLASMIDS_HOST (
        final_metadata.map { r -> record(meta: r.plasmid_meta, table: r.plasmids_final_metadata) }
    )

    // Plot Sankey for hosts
    SANKEY_VIRUSES_HOST (
        final_metadata.map { r -> record(meta: r.viral_meta, table: r.viruses_final_metadata) }
    )

    //
    // Collate and save software versions
    //
    // Processes send (process, tool, version) to the `versions` topic. Rebuild
    // the per-process YAML that versions.yml files used to hold.
    def versions = channel.topic('versions') as Channel<Tuple<String,String,String>>
    collated_versions = versions
        .unique()
        .map { process, tool, version -> tuple(process.tokenize(':')[-1], "    ${tool}: ${version}") }
        .groupBy()
        .map { process, tool_versions -> processVersionsFromYAML("\"${process}\":\n" + tool_versions.toUnique().toSorted().join('\n')) }
        .unique()
        .mix(channel.of(workflowVersionToYAML()))
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'metaviraverse_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        channel.fromPath(params.multiqc_config, checkIfExists: true) :
        channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = channel.of(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_workflow_summary
        .collectFile(name: 'workflow_summary_mqc.yaml')
        .mix(collated_versions)

    multiqc = MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.collect(),
        ch_multiqc_custom_config.collect(),
        ch_multiqc_logo.collect(),
        [],
        []
    )
    multiqc_report = multiqc.report.collect() // channel: /path/to/multiqc_report.html

    emit:
    multiqc_report = multiqc_report

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
