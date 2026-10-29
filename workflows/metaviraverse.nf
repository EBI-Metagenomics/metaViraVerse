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

include { PREPROCESSING ; PreprocessingParams ; Preprocessed } from '../subworkflows/local/preprocessing/main'
include { PROCESS_VIRAL_SEQUENCES ; ViralParams ; ViralCatalogue } from '../subworkflows/local/process_viral_sequences'
include { PROCESS_PLASMIDS ; PlasmidParams ; PlasmidCatalogue } from '../subworkflows/local/process_plasmids'
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

record MetaviraverseParams {
    catalogues_metadata: Path
    multiqc_config: Path?
    multiqc_logo: String?
}

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

record CatalogueSummary {
    catalogue_json: Path            // statistics of the whole catalogue
    viral_meta: Map
    plasmid_meta: Map
    viruses_metadata: Path
    viruses_cluster_stats: Path
    plasmids_metadata: Path
    viruses_source_sankey: Path?    // taxonomy of the source genomes
    plasmids_source_sankey: Path?
}

workflow METAVIRAVERSE {

    take:
    samples: Channel<Sample>                       // samplesheet read in from --input
    third_party_input: Channel<ThirdPartyInput>    // samplesheet read in from --third_party_input (may be empty)
    opts: MetaviraverseParams
    preprocessing_opts: PreprocessingParams
    viral_opts: ViralParams
    plasmid_opts: PlasmidParams

    main:

    //
    // ----------- Preprocess third party data (coming not from MGnify)
    //
    // (runs no tasks when there is no third-party input)
    third_party_data = THIRD_PARTY_DATA(
        third_party_input
    )
    third_party = third_party_data.third_party_data
    ch_invalid_fastas = third_party_data.invalid_fastas

    //
    // Separate viral sequences and plasmids
    //
    preprocessed = PREPROCESSING (
       samples,
       third_party,
       preprocessing_opts
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
       preprocessed.map { r -> r.map_file },
       viral_opts
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
       preprocessed.map { r -> r.filtered_gff },  // concatenated GFF built from all inputs
       plasmid_opts
    )

    // one record with both catalogues and the preprocessing results
    catalogues = viral_catalogue
        .combine(plasmid_catalogue)
        .map { viral, plasmid -> record(viral: viral, plasmid: plasmid) }
        .combine(preprocessed: preprocessed)

    //
    // Collect stats for whole catalogue into JSON
    //
    catalogue_stats = COLLECT_CATALOGUE_STATS (
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
                vitap_best: r.viral.taxonomy.vitap_best,
                viphogs_taxonomy: r.viral.taxonomy.viphogs,
                genomad: r.viral.taxonomy.genomad,
                map_file: r.preprocessed.map_file,
                iphop_genome: r.viral.hosts.iphop_genome,
                iphop_genus: r.viral.hosts.iphop_genus,
                spacepharer: r.viral.hosts.spacepharer
            )
        },
        opts.catalogues_metadata
    )


    // Plot Sankey for hosts
    plasmids_sankey = SANKEY_PLASMIDS_HOST (
        final_metadata.map { r -> record(meta: r.plasmid_meta, table: r.plasmids_final_metadata) }
    )

    // Plot Sankey for hosts
    viruses_sankey = SANKEY_VIRUSES_HOST (
        final_metadata.map { r -> record(meta: r.viral_meta, table: r.viruses_final_metadata) }
    )

    // one record for the whole catalogue (each channel has a single item)
    summary = final_metadata
        .combine(catalogue_stats)
        .map { m, s ->
            record(
                catalogue_json: s.catalogue_json,
                viral_meta: m.viral_meta,
                plasmid_meta: m.plasmid_meta,
                viruses_metadata: m.viruses_final_metadata,
                viruses_cluster_stats: m.viruses_reps_stats,
                plasmids_metadata: m.plasmids_final_metadata
            )
        }
        .combine(viruses_sankey)
        .map { r, v -> r + record(viruses_source_sankey: v.sankey_html) }
        .combine(plasmids_sankey)
        .map { r, p -> r + record(plasmids_source_sankey: p.sankey_html) }

    //
    // Collate and save software versions
    //
    // Processes send (process, tool, version) to the `versions` topic. Rebuild
    // the per-process YAML that versions.yml files used to hold.
    def version_topic = channel.topic('versions') as Channel<Tuple<String,String,String>>
    collated_versions = version_topic
        .unique()
        .map { process, tool, version -> tuple(process.tokenize(':')[-1], "    ${tool}: ${version}") }
        .groupBy()
        .map { process, tool_versions -> processVersionsFromYAML("\"${process}\":\n" + tool_versions.toUnique().toSorted().join('\n')) }
        .unique()
        .mix(channel.of(workflowVersionToYAML()))
        .collectFile(
            name:  'metaviraverse_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = opts.multiqc_config ?
        channel.of(opts.multiqc_config) :
        channel.empty()
    ch_multiqc_logo          = opts.multiqc_logo ?
        channel.fromPath(opts.multiqc_logo, checkIfExists: true) :
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
    multiqc_files = multiqc.report.mix(multiqc.data, multiqc.plots)   // report, data and (optional) plots

    emit:
    multiqc_report = multiqc_report
    preprocessed: Value<Preprocessed> = preprocessed
    viruses: Channel<ViralCatalogue> = viral_catalogue
    plasmids: Channel<PlasmidCatalogue> = plasmid_catalogue
    summary: Channel<CatalogueSummary> = summary
    multiqc: Channel<Path> = multiqc_files
    versions: Channel<Path> = collated_versions
    invalid_fastas: Channel<Path> = ch_invalid_fastas   // third-party FASTAs that failed validation

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
