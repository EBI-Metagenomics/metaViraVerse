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

include { METAVIRAVERSE ; MetaviraverseParams ; CatalogueSummary } from './workflows/metaviraverse'
include { PreprocessingParams ; Preprocessed } from './subworkflows/local/preprocessing/main'
include { ViralParams ; ViralCatalogue } from './subworkflows/local/process_viral_sequences'
include { PlasmidParams ; PlasmidCatalogue } from './subworkflows/local/process_plasmids'
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
    input: Channel<SamplesheetRow>?                 // MGnify samplesheet (CSV, JSON or YAML)
    third_party_input: Channel<ThirdPartyRow>?      // third-party samplesheet
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
    viral_clustering: String = 'default'   // 'default' or 'vjdb' (pipelines/vjdb-clustering)

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

    // Publishing
    save_intermediates: Boolean = false

    // Boilerplate
    monochrome_logs: Boolean = false
    version: Boolean = false
    validate_params: Boolean = true
}

// one row of --input
record SamplesheetRow {
    id: String
    fna: Path
    gff: Path
    faa: Path
    source: String?     // genome | metagenome
    biome: String?
}

// one row of --third_party_input
record ThirdPartyRow {
    id: String
    fasta: Path
    type: String        // virus | plasmid | prophage
    source: String?     // genome | metagenome
    biome: String?
    study_accession: String?
    sample_accession: String?
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
    samplesheet: Channel<SamplesheetRow>
    samplesheet_third_party: Channel<ThirdPartyRow>
    opts: MetaviraverseParams
    preprocessing_opts: PreprocessingParams
    viral_opts: ViralParams
    plasmid_opts: PlasmidParams

    main:

    //
    // WORKFLOW: Run pipeline
    //
    // Rows become [meta, files] records. An empty optional column is passed on as '[]', the
    // string the samplesheetToList value (an empty list) used to become in the rename map.
    metaviraverse = METAVIRAVERSE (
        samplesheet.map { s ->
            checkSamplesheetRow(s)
            record(meta: [id: s.id, source: s.source ?: '[]', biome: s.biome ?: '[]'], fna: s.fna, gff: s.gff, faa: s.faa)
        },
        samplesheet_third_party.map { s ->
            checkThirdPartyRow(s)
            def meta = [
                id: s.id, type: s.type, source: s.source ?: '[]', biome: s.biome ?: '[]',
                study_accession: s.study_accession ?: '[]', sample_accession: s.sample_accession ?: '[]'
            ]
            record(meta: meta, fasta: s.fasta)
        },
        opts,
        preprocessing_opts,
        viral_opts,
        plasmid_opts
    )
    multiqc_report = metaviraverse.multiqc_report
    preprocessed = metaviraverse.preprocessed
    viruses = metaviraverse.viruses
    plasmids = metaviraverse.plasmids
    summary = metaviraverse.summary
    multiqc = metaviraverse.multiqc
    versions = metaviraverse.versions
    invalid_fastas = metaviraverse.invalid_fastas

    emit:
    multiqc_report = multiqc_report // channel: /path/to/multiqc_report.html
    preprocessed: Value<Preprocessed> = preprocessed
    viruses: Channel<ViralCatalogue> = viruses
    plasmids: Channel<PlasmidCatalogue> = plasmids
    summary: Channel<CatalogueSummary> = summary
    multiqc: Channel<Path> = multiqc
    versions: Channel<Path> = versions
    invalid_fastas: Channel<Path> = invalid_fastas
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
        params.outdir
    )

    //
    // WORKFLOW: Run main workflow
    //
    metaviraverse = EBIMETAGENOMICS_METAVIRAVERSE (
        params.input ?: channel.empty(),
        params.third_party_input ?: channel.empty(),
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
        metaviraverse.multiqc_report
    )

    publish:
    preprocessed = metaviraverse.preprocessed
    preprocessed_intermediates = metaviraverse.preprocessed
    viruses = metaviraverse.viruses
    plasmids = metaviraverse.plasmids
    plasmids_intermediates = metaviraverse.plasmids
    summary = metaviraverse.summary
    multiqc = metaviraverse.multiqc
    versions = metaviraverse.versions
    invalid_fastas = metaviraverse.invalid_fastas
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Samplesheet checks beyond the record types (formerly assets/schema_input.json)
def checkSamplesheetRow(s: SamplesheetRow) {
    if( !(s.id ==~ /\S+/) )
        error("Invalid --input row '${s.id}': id must not contain spaces")
    if( s.source && !(s.source in ['genome', 'metagenome']) )
        error("Invalid --input row '${s.id}': source must be 'genome' or 'metagenome', got '${s.source}'")
    if( s.biome && !(s.biome ==~ /\S+/) )
        error("Invalid --input row '${s.id}': biome must not contain spaces")
}

// (formerly assets/schema_third_party_input.json)
def checkThirdPartyRow(s: ThirdPartyRow) {
    if( !(s.id ==~ /\S+/) )
        error("Invalid --third_party_input row '${s.id}': id must not contain spaces")
    if( !(s.type in ['virus', 'plasmid', 'prophage']) )
        error("Invalid --third_party_input row '${s.id}': type must be 'virus', 'plasmid' or 'prophage', got '${s.type}'")
    if( s.source && !(s.source in ['genome', 'metagenome']) )
        error("Invalid --third_party_input row '${s.id}': source must be 'genome' or 'metagenome', got '${s.source}'")
    if( s.biome && !(s.biome ==~ /\S+/) )
        error("Invalid --third_party_input row '${s.id}': biome must not contain spaces")
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW OUTPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Published into `outputDir` (--outdir); see docs/output.md for the layout.
    `*_intermediates` outputs are published only with --save_intermediates.
    Null fields (skipped steps, missing optional files) are not published.
*/

output {
    // rename map, catalogue metadata tables and the chosen sequences (before clustering)
    preprocessed: Preprocessed {
        path { r ->
            r.map_file >> "additional_data/rename_mapfile.tsv"
            r.excluded_metadata >> "additional_data/catalogue_excluded_records_metadata.tsv"
            r.metadata >> "additional_data/catalogue_filtered_metadata.tsv"
            r.virus_gz.fna >> "viruses/viral_sequences.fna.gz"
            r.virus_gz.gff >> "viruses/viral_sequences.gff.gz"
            r.virus_gz.faa >> "viruses/viral_sequences.faa.gz"
            r.prophage_gz.fna >> "viruses/prophages.fna.gz"
            r.prophage_gz.gff >> "viruses/prophages.gff.gz"
            r.prophage_gz.faa >> "viruses/prophages.faa.gz"
            r.plasmid_gz.fna >> "plasmids/plasmids.fna.gz"
            r.plasmid_gz.gff >> "plasmids/plasmids.gff.gz"
            r.plasmid_gz.faa >> "plasmids/plasmids.faa.gz"
        }
    }

    preprocessed_intermediates: Preprocessed {
        enabled params.save_intermediates
        path { r ->
            r.fna_renamed >> "additional_data/rename/"
            r.gff_renamed >> "additional_data/rename/"
            r.checkv_report >> "additional_data/viruses_checkv_report.tsv"
            r.rrna_gff >> "additional_data/barrnap/"
            r.filtered_metadata >> "additional_data/choose_sequences/"
            r.filtered_fna >> "additional_data/choose_sequences/"
            r.filtered_gff >> "additional_data/choose_sequences/"
            r.filtered_faa >> "additional_data/choose_sequences/"
            r.virus_fna >> "additional_data/choose_sequences/"
            r.virus_gff >> "additional_data/choose_sequences/"
            r.virus_faa >> "additional_data/choose_sequences/"
            r.prophage_fna >> "additional_data/choose_sequences/"
            r.prophage_gff >> "additional_data/choose_sequences/"
            r.prophage_faa >> "additional_data/choose_sequences/"
            r.plasmid_fna >> "additional_data/choose_sequences/"
            r.plasmid_gff >> "additional_data/choose_sequences/"
            r.plasmid_faa >> "additional_data/choose_sequences/"
        }
    }

    // viral cluster representatives and their annotations
    viruses: Channel<ViralCatalogue> {
        path { r ->
            // CheckV ANI clustering (--cluster_vclust false)
            r.ani_tsv >> "${r.meta.id}/clustering/"
            r.ani_clusters_tsv >> "${r.meta.id}/clustering/"
            // vjdb-clustering membership table (--viral_clustering vjdb)
            r.vjdb_clusters_csv >> "${r.meta.id}/clustering/"

            r.final_gff >> "${r.meta.id}/cluster_representatives/viruses_reps_annotated.gff.gz"
            r.index.fna_gz >> "${r.meta.id}/cluster_representatives/"
            r.index.fna_fai >> "${r.meta.id}/cluster_representatives/"
            r.index.fna_gzi >> "${r.meta.id}/cluster_representatives/"
            r.index.faa_gz >> "${r.meta.id}/cluster_representatives/"
            r.index.faa_fai >> "${r.meta.id}/cluster_representatives/"
            r.index.faa_gzi >> "${r.meta.id}/cluster_representatives/"
            r.index.gff_gz >> "${r.meta.id}/cluster_representatives/"
            r.index.gff_tbi >> "${r.meta.id}/cluster_representatives/"
            r.index.bacphlip_gz >> "${r.meta.id}/cluster_representatives/functional_annotation/lifestyle/"
            r.index.bacphlip_tbi >> "${r.meta.id}/cluster_representatives/functional_annotation/lifestyle/"

            // taxonomy: geNomad
            r.taxonomy.genomad >> "${r.meta.id}/cluster_representatives/taxonomy/genomad/"
            r.taxonomy.genomad_counts >> "${r.meta.id}/cluster_representatives/taxonomy/genomad/"
            r.taxonomy.genomad_plots?.krona >> "${r.meta.id}/cluster_representatives/taxonomy/genomad/"
            r.taxonomy.genomad_plots?.sankey >> "${r.meta.id}/cluster_representatives/taxonomy/genomad/"
            r.taxonomy.genomad_plots?.itol_tree >> "${r.meta.id}/cluster_representatives/taxonomy/genomad/itol_genomad/"
            r.taxonomy.genomad_plots?.itol_labels >> "${r.meta.id}/cluster_representatives/taxonomy/genomad/itol_genomad/"
            r.taxonomy.genomad_plots?.itol_counts >> "${r.meta.id}/cluster_representatives/taxonomy/genomad/itol_genomad/"
            r.taxonomy.genomad_plots?.itol_types >> "${r.meta.id}/cluster_representatives/taxonomy/genomad/itol_genomad/"
            r.taxonomy.genomad_plots?.itol_biomes >> "${r.meta.id}/cluster_representatives/taxonomy/genomad/itol_genomad/"
            // taxonomy: ViPhOGs
            r.taxonomy.viphogs >> "${r.meta.id}/cluster_representatives/taxonomy/viphogs/"
            r.taxonomy.viphogs_modified >> "${r.meta.id}/cluster_representatives/taxonomy/viphogs/"
            r.taxonomy.viphogs_informative >> "${r.meta.id}/cluster_representatives/taxonomy/viphogs/"
            r.taxonomy.viphogs_annotation >> "${r.meta.id}/cluster_representatives/taxonomy/viphogs/"
            r.taxonomy.viphogs_counts >> "${r.meta.id}/cluster_representatives/taxonomy/viphogs/"
            r.taxonomy.viphogs_plots?.krona >> "${r.meta.id}/cluster_representatives/taxonomy/viphogs/"
            r.taxonomy.viphogs_plots?.sankey >> "${r.meta.id}/cluster_representatives/taxonomy/viphogs/"
            r.taxonomy.viphogs_plots?.itol_tree >> "${r.meta.id}/cluster_representatives/taxonomy/viphogs/itol_viphogs/"
            r.taxonomy.viphogs_plots?.itol_labels >> "${r.meta.id}/cluster_representatives/taxonomy/viphogs/itol_viphogs/"
            r.taxonomy.viphogs_plots?.itol_counts >> "${r.meta.id}/cluster_representatives/taxonomy/viphogs/itol_viphogs/"
            r.taxonomy.viphogs_plots?.itol_types >> "${r.meta.id}/cluster_representatives/taxonomy/viphogs/itol_viphogs/"
            r.taxonomy.viphogs_plots?.itol_biomes >> "${r.meta.id}/cluster_representatives/taxonomy/viphogs/itol_viphogs/"
            // taxonomy: VITAP
            r.taxonomy.vitap_best >> "${r.meta.id}/cluster_representatives/taxonomy/vitap/"
            r.taxonomy.vitap_counts >> "${r.meta.id}/cluster_representatives/taxonomy/vitap/"
            r.taxonomy.vitap_plots?.krona >> "${r.meta.id}/cluster_representatives/taxonomy/vitap/"
            r.taxonomy.vitap_plots?.sankey >> "${r.meta.id}/cluster_representatives/taxonomy/vitap/"
            r.taxonomy.vitap_plots?.itol_tree >> "${r.meta.id}/cluster_representatives/taxonomy/vitap/itol_vitap/"
            r.taxonomy.vitap_plots?.itol_labels >> "${r.meta.id}/cluster_representatives/taxonomy/vitap/itol_vitap/"
            r.taxonomy.vitap_plots?.itol_counts >> "${r.meta.id}/cluster_representatives/taxonomy/vitap/itol_vitap/"
            r.taxonomy.vitap_plots?.itol_types >> "${r.meta.id}/cluster_representatives/taxonomy/vitap/itol_vitap/"
            r.taxonomy.vitap_plots?.itol_biomes >> "${r.meta.id}/cluster_representatives/taxonomy/vitap/itol_vitap/"

            // host prediction: iPHoP
            r.hosts.iphop_genome >> "${r.meta.id}/cluster_representatives/host_detection/iphop/${r.hosts.iphop_genome?.name}.csv"
            r.hosts.iphop_genus >> "${r.meta.id}/cluster_representatives/host_detection/iphop/${r.hosts.iphop_genus?.name}.csv"
            r.hosts.iphop_genome_sankey >> "${r.meta.id}/cluster_representatives/host_detection/iphop/"
            r.hosts.iphop_genome_reduced_sankey >> "${r.meta.id}/cluster_representatives/host_detection/iphop/"
            r.hosts.iphop_genus_sankey >> "${r.meta.id}/cluster_representatives/host_detection/iphop/"
            r.hosts.iphop_genus_reduced_sankey >> "${r.meta.id}/cluster_representatives/host_detection/iphop/"
            // host prediction: SpacePHARER (--predict_host_from_custom_spacers)
            r.hosts.spacepharer_predictions >> "${r.meta.id}/cluster_representatives/host_detection/spacepharer/"
            r.hosts.spacepharer >> "${r.meta.id}/cluster_representatives/host_detection/spacepharer/"
            r.hosts.spacepharer_lineage_json >> "${r.meta.id}/cluster_representatives/host_detection/spacepharer/"
            r.hosts.spacepharer_sankey >> "${r.meta.id}/cluster_representatives/host_detection/spacepharer/"
            r.hosts.spacepharer_reduced_sankey >> "${r.meta.id}/cluster_representatives/host_detection/spacepharer/"

            // proteins: AMR (one report per tool), HMMER, phams
            r.proteins.amrfinderplus_tsv >> "${r.meta.id}/cluster_representatives/functional_annotation/amr/${r.meta.id}_amrfinderplus.tsv"
            r.proteins.deeparg_tsv >> "${r.meta.id}/cluster_representatives/functional_annotation/amr/${r.meta.id}_deeparg.tsv"
            r.proteins.rgi_tsv >> "${r.meta.id}/cluster_representatives/functional_annotation/amr/${r.meta.id}_rgi.tsv"
            r.proteins.amr_gff >> "${r.meta.id}/cluster_representatives/functional_annotation/amr/"
            r.proteins.hmmer_tables >> "${r.meta.id}/cluster_representatives/functional_annotation/hmmer/"
            r.proteins.hmmer_summaries >> "${r.meta.id}/cluster_representatives/functional_annotation/hmmer/"
            r.proteins.phams_tsv >> "${r.meta.id}/cluster_representatives/protein_clusters_phammseqs/"
            r.proteins.phams_plots?.size_distribution_tsv >> "${r.meta.id}/cluster_representatives/protein_clusters_phammseqs/"
            r.proteins.phams_plots?.size_distribution_png >> "${r.meta.id}/cluster_representatives/protein_clusters_phammseqs/"
            r.proteins.phams_plots?.rank_abundance >> "${r.meta.id}/cluster_representatives/protein_clusters_phammseqs/"
            r.proteins.phams_plots?.top_functions_tsv >> "${r.meta.id}/cluster_representatives/protein_clusters_phammseqs/"
            r.proteins.phams_plots?.top_functions_png >> "${r.meta.id}/cluster_representatives/protein_clusters_phammseqs/"
            r.proteins.phams_plots?.purity_tsv >> "${r.meta.id}/cluster_representatives/protein_clusters_phammseqs/"
            r.proteins.phams_plots?.purity_png >> "${r.meta.id}/cluster_representatives/protein_clusters_phammseqs/"
            r.proteins.phams_plots?.stats >> "${r.meta.id}/cluster_representatives/protein_clusters_phammseqs/"
        }
    }

    // plasmid cluster representatives and their annotations
    plasmids: Channel<PlasmidCatalogue> {
        path { r ->
            // CheckV ANI clustering (--cluster_vclust false)
            r.ani_tsv >> "${r.meta.id}/clustering/"
            r.ani_clusters_tsv >> "${r.meta.id}/clustering/"

            r.annotated_gff >> "${r.meta.id}/cluster_representatives/plasmids_annotated.gff.gz"
            r.index.fna_gz >> "${r.meta.id}/cluster_representatives/"
            r.index.fna_fai >> "${r.meta.id}/cluster_representatives/"
            r.index.fna_gzi >> "${r.meta.id}/cluster_representatives/"
            r.index.faa_gz >> "${r.meta.id}/cluster_representatives/"
            r.index.faa_fai >> "${r.meta.id}/cluster_representatives/"
            r.index.faa_gzi >> "${r.meta.id}/cluster_representatives/"
            r.index.gff_gz >> "${r.meta.id}/cluster_representatives/"
            r.index.gff_tbi >> "${r.meta.id}/cluster_representatives/"

            // mobility: plaSquid + MOB-suite
            r.mobility_stats >> "${r.meta.id}/cluster_representatives/functional_annotation/"
            r.plasquid.report >> "${r.meta.id}/cluster_representatives/functional_annotation/plasquid/final_report_per_contig.tsv"
            r.plasquid.protein_report >> "${r.meta.id}/cluster_representatives/functional_annotation/plasquid/final_report_per_protein.tsv"
            r.plasquid.mobility_classification >> "${r.meta.id}/cluster_representatives/functional_annotation/plasquid/"
            r.plasquid.mobility_stats >> "${r.meta.id}/cluster_representatives/functional_annotation/plasquid/mobility_stats_base.json"
            r.mobsuite.biomarker_report >> "${r.meta.id}/cluster_representatives/functional_annotation/mobsuite-typer/"
            r.mobsuite.mge_report >> "${r.meta.id}/cluster_representatives/functional_annotation/mobsuite-typer/"
            r.mobsuite.mob_report >> "${r.meta.id}/cluster_representatives/functional_annotation/mobsuite-typer/"

            // proteins: AMR (one report per tool), MMseqs2 clusters
            r.proteins.amrfinderplus_tsv >> "${r.meta.id}/cluster_representatives/functional_annotation/amr/${r.meta.id}_amrfinderplus.tsv"
            r.proteins.deeparg_tsv >> "${r.meta.id}/cluster_representatives/functional_annotation/amr/${r.meta.id}_deeparg.tsv"
            r.proteins.rgi_tsv >> "${r.meta.id}/cluster_representatives/functional_annotation/amr/${r.meta.id}_rgi.tsv"
            r.proteins.amr_gff >> "${r.meta.id}/cluster_representatives/functional_annotation/amr/"
            r.proteins.mmseqs_tsv >> "${r.meta.id}/cluster_representatives/protein_clusters_mmseqs/"
            r.proteins.mmseqs_plots?.size_distribution_tsv >> "${r.meta.id}/cluster_representatives/protein_clusters_mmseqs/"
            r.proteins.mmseqs_plots?.size_distribution_png >> "${r.meta.id}/cluster_representatives/protein_clusters_mmseqs/"
            r.proteins.mmseqs_plots?.rank_abundance >> "${r.meta.id}/cluster_representatives/protein_clusters_mmseqs/"
            r.proteins.mmseqs_plots?.top_functions_tsv >> "${r.meta.id}/cluster_representatives/protein_clusters_mmseqs/"
            r.proteins.mmseqs_plots?.top_functions_png >> "${r.meta.id}/cluster_representatives/protein_clusters_mmseqs/"
            r.proteins.mmseqs_plots?.purity_tsv >> "${r.meta.id}/cluster_representatives/protein_clusters_mmseqs/"
            r.proteins.mmseqs_plots?.purity_png >> "${r.meta.id}/cluster_representatives/protein_clusters_mmseqs/"
            r.proteins.mmseqs_plots?.stats >> "${r.meta.id}/cluster_representatives/protein_clusters_mmseqs/"
        }
    }

    // plaSquid evidence tables
    plasmids_intermediates: Channel<PlasmidCatalogue> {
        enabled params.save_intermediates
        path { r ->
            r.plasquid.mob_candidates >> "${r.meta.id}/cluster_representatives/functional_annotation/plasquid/evidence/mobilisation_candidates.tsv"
            r.plasquid.rep_candidates >> "${r.meta.id}/cluster_representatives/functional_annotation/plasquid/evidence/replication_initiator_proteins.tsv"
            r.plasquid.rna_candidates >> "${r.meta.id}/cluster_representatives/functional_annotation/plasquid/evidence/"
            r.plasquid.inc_candidates >> "${r.meta.id}/cluster_representatives/functional_annotation/plasquid/evidence/incompatibility_candidates.tsv"
        }
    }

    // catalogue statistics, final metadata tables and source taxonomy plots
    summary: Channel<CatalogueSummary> {
        path { r ->
            r.catalogue_json >> "./"
            r.viruses_metadata >> "viruses/"
            r.viruses_cluster_stats >> "viruses/"
            r.plasmids_metadata >> "plasmids/"
            r.viruses_source_sankey >> "${r.viral_meta.id}/"
            r.plasmids_source_sankey >> "${r.plasmid_meta.id}/"
        }
    }

    multiqc: Channel<Path> {
        enabled params.save_intermediates
        path 'multiqc'
    }

    versions: Channel<Path> {
        path 'pipeline_info'
    }

    // third-party FASTAs that failed validation
    invalid_fastas: Channel<Path> {
        path '.'
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
