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
include { GENERATE_README                       } from '../modules/local/generate_readme'
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
    analyses_metadata: Path?   // fills source study/sample of sequences named <ERZ accession>|...
    multiqc_config: Path?
    multiqc_logo: String?
    // title section of README.md
    catalogue_name: String?
    catalogue_version: String?
    catalogue_website_url: String?
    source_catalogues: String?  // comma-separated URLs
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
    viruses_source_sankey_table: Path?    // links the source Sankey plots are drawn from
    plasmids_source_sankey_table: Path?
    readme: Path?                   // README.md of the catalogue
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
                map_file: r.preprocessed.map_file,
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
        opts.catalogues_metadata,
        opts.analyses_metadata
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
    catalogue_summary = final_metadata
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
        .map { r, v -> r + record(viruses_source_sankey: v.sankey_html, viruses_source_sankey_table: v.sankey_tsv) }
        .combine(plasmids_sankey)
        .map { r, p -> r + record(plasmids_source_sankey: p.sankey_html, plasmids_source_sankey_table: p.sankey_tsv) }

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
    // MODULE: README.md of the catalogue
    //
    // Options that decide which tools ran, for the tool list and methods text of the README
    run_options = third_party_input.count().map { n ->
        groovy.json.JsonOutput.toJson([
            third_party_input: n > 0,
            skip_rrna_detection: preprocessing_opts.skip_rrna_detection,
            skip_checkv: preprocessing_opts.skip_checkv,
            rename_accession: preprocessing_opts.rename_accession,
            viral_clustering: viral_opts.viral_clustering,
            cluster_vclust: viral_opts.cluster_vclust,
            skip_vitap: viral_opts.skip_vitap,
            skip_genomad: viral_opts.skip_genomad,
            skip_viphogs: viral_opts.skip_viphogs,
            skip_iphop: viral_opts.skip_iphop,
            predict_host_from_custom_spacers: viral_opts.predict_host_from_custom_spacers,
            skip_phammseqs: viral_opts.skip_phammseqs,
            skip_amrfinderplus_for_viruses: viral_opts.skip_amrfinderplus_for_viruses,
            skip_deeparg_for_viruses: viral_opts.skip_deeparg_for_viruses,
            skip_rgi_for_viruses: viral_opts.skip_rgi_for_viruses,
            skip_mmseqs: plasmid_opts.skip_mmseqs,
            mmseqs_id_threshold: plasmid_opts.mmseqs_id_threshold,
            mmseqs_cov_threshold: plasmid_opts.mmseqs_cov_threshold,
            skip_amrfinderplus_for_plasmids: plasmid_opts.skip_amrfinderplus_for_plasmids,
            skip_deeparg_for_plasmids: plasmid_opts.skip_deeparg_for_plasmids,
            skip_rgi_for_plasmids: plasmid_opts.skip_rgi_for_plasmids
        ])
    }

    readme = GENERATE_README (
        catalogue_summary
            .combine(collated_versions)
            .map { r, versions ->
                record(
                    catalogue_json: r.catalogue_json,
                    viruses_metadata: r.viruses_metadata,
                    plasmids_metadata: r.plasmids_metadata,
                    versions_yml: versions
                )
            },
        run_options,
        opts.catalogue_name,
        opts.catalogue_version,
        opts.catalogue_website_url,
        opts.source_catalogues
    )

    summary = catalogue_summary
        .combine(readme)
        .map { r, m -> r + record(readme: m.readme) }

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
