/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap                      } from 'plugin/nf-schema'
include { paramsSummaryMultiqc                  } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                } from '../subworkflows/local/utils_nfcore_metaviraverse_pipeline'

include { PREPROCESSING                         } from '../subworkflows/local/preprocessing'
include { PROCESS_VIRAL_SEQUENCES               } from '../subworkflows/local/process_viral_sequences'
include { PROCESS_PLASMIDS                      } from '../subworkflows/local/process_plasmids'

include { COLLECT_CATALOGUE_STATS               } from '../modules/local/collect_catalogue_stats'

include { MULTIQC                               } from '../modules/nf-core/multiqc'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow METAVIRAVERSE {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    //
    // Separate viral sequences and plasmids
    //
    PREPROCESSING (
       ch_samplesheet
    )
    ch_versions = ch_versions.mix(PREPROCESSING.out.versions)

    // Process viral_sequences and prophages together as "viruses"

    viruses = PREPROCESSING.out.viral_sequences
        .map{ _meta, seqs -> seqs }
        .combine(PREPROCESSING.out.prophages
        .map{ _meta, seqs -> seqs })
        .flatMap { tuple -> tuple }
        .collectFile(name: "viruses.fasta")
        .map{ seqs -> [[id: 'viruses'], seqs]}

    // publish viruses
    viruses.subscribe { meta, seqs ->
        def outDir = file("${params.outdir}/${meta.id}")
        outDir.mkdirs()  // Create directory if it doesn't exist

        def outPath = file("${outDir}/viruses.fasta.gz")
        outPath.withOutputStream { out ->
            new java.util.zip.GZIPOutputStream(out).withWriter { writer ->
                writer << seqs.text
            }
        }
    }

    //
    // Process viruses
    //
    PROCESS_VIRAL_SEQUENCES (
       viruses,
       PREPROCESSING.out.combined_gff,   // concatenated GFF built from all inputs
       PREPROCESSING.out.combined_faa,   // concatenated FAA built from all inputs
       PREPROCESSING.out.metadata,
       []
    )
    ch_versions = ch_versions.mix(PROCESS_VIRAL_SEQUENCES.out.versions)

    // Process plasmids filtered from input

    plasmids = PREPROCESSING.out.plasmids
       .map{ _meta, seqs -> seqs }
       .collectFile(name: "plasmids.fasta")
       .map{ seqs -> [[id: 'plasmids'], seqs]}

    // publish plasmids
    plasmids.subscribe { meta, seqs ->
        def outDir = file("${params.outdir}/${meta.id}")
        outDir.mkdirs()  // Create directory if it doesn't exist

        def outPath = file("${outDir}/plasmids.fasta.gz")
        outPath.withOutputStream { out ->
            new java.util.zip.GZIPOutputStream(out).withWriter { writer ->
                writer << seqs.text
            }
        }
    }

    //
    // Process plasmids
    //
    PROCESS_PLASMIDS (
       plasmids,
       PREPROCESSING.out.combined_faa,
       PREPROCESSING.out.combined_gff,   // concatenated GFF built from all inputs
       []
    )
    ch_versions = ch_versions.mix(PROCESS_PLASMIDS.out.versions)


    //
    // Collect stats for whole catalogue into JSON
    //
    COLLECT_CATALOGUE_STATS (
        PREPROCESSING.out.viral_sequences,
        PREPROCESSING.out.prophages,
        PREPROCESSING.out.plasmids,
        PREPROCESSING.out.metadata,
        PROCESS_VIRAL_SEQUENCES.out.reps_tsv,
        PROCESS_PLASMIDS.out.reps_tsv,
        PROCESS_VIRAL_SEQUENCES.out.reps_proteins,
        PROCESS_PLASMIDS.out.reps_proteins
    )

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'metaviraverse_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

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
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
