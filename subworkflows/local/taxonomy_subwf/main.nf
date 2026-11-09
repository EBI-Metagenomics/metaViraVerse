nextflow.enable.types = true

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
include { TAXONOMY_VISUALISATION as VIS_VIPHOGS ; TaxonomyPlots } from '../taxonomy_visualisation/main'
include { TAXONOMY_VISUALISATION as VIS_VITAP          } from '../taxonomy_visualisation/main'
include { TAXONOMY_VISUALISATION as VIS_GENOMAD        } from '../taxonomy_visualisation/main'



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

record TaxonomyParams {
    skip_vitap: Boolean
    skip_genomad: Boolean
    skip_viphogs: Boolean
    vitap_db: Path?
    genomad_db: Path?
    // ViPhOGs (see ViphogsParams)
    viphog_db: Path?
    additional_viphog_model_data: Path?
    ncbi_db: Path?
    viphogs_factor_file: Path?
    viphogs_evalue: Float
    protein_annotation_fasta_chunksize: Integer
}

record TaxonomyInput {
    meta: Map
    faa: Path        // representative proteins, uncompressed
    reps_gff: Path
}

record FnaChunk {
    meta: Map
    fasta: Path
}

record VitapBest {
    meta: Map
    vitap_best: Path
    vitap_counts: Path
    vitap_plots: TaxonomyPlots
}

record GenomadTaxonomy {
    meta: Map
    genomad: Path
    genomad_counts: Path
    genomad_plots: TaxonomyPlots
}

record ViphogsTaxonomy {
    meta: Map
    viphogs: Path
    viphogs_modified: Path
    viphogs_informative: Path
    viphogs_annotation: Path
    viphogs_counts: Path
    viphogs_plots: TaxonomyPlots
}

// Fields of a skipped tool are null
record TaxonomyResult {
    meta: Map
    vitap_best: Path?
    vitap_counts: Path?            // taxonomy counts joined with metadata
    vitap_plots: TaxonomyPlots?
    genomad: Path?
    genomad_counts: Path?
    genomad_plots: TaxonomyPlots?
    viphogs: Path?                 // ViPhOGs taxonomy assignment
    viphogs_modified: Path?
    viphogs_informative: Path?
    viphogs_annotation: Path?
    viphogs_counts: Path?
    viphogs_plots: TaxonomyPlots?
}

// taxonomy tables and plots carry the tool in their meta: drop it to join with the input
def withoutTool(meta: Map) -> Map {
    return meta.subMap(meta.keySet().findAll { k -> k != 'tool' })
}

workflow TAXONOMY_ASSIGNMENT {

    take:
    reps: Channel<TaxonomyInput>
    fna_chunks: Channel<FnaChunk>
    combined_metadata: Value<Path>
    opts: TaxonomyParams

    main:

    def vitap_best: Channel<VitapBest> = channel.empty()
    def genomad_taxonomy: Channel<GenomadTaxonomy> = channel.empty()
    def viphogs_taxonomy: Channel<ViphogsTaxonomy> = channel.empty()

    def ch_fna_chunks = fna_chunks.map { r -> tuple(r.meta, r.fasta) }

    //
    // Taxonomy genoMAD
    //
    if ( !opts.skip_genomad ) {
        genomad = GENOMAD_ENDTOEND (
            ch_fna_chunks,
            opts.genomad_db
        )

        genomad_concat = CONCATENATE_GENOMAD (
            genomad.virus_summary.groupTuple(),
            'tsv',
            'tsv'
        )

        tax_genomad = TAX_GENOMAD (
            genomad_concat.csv.map { meta, table -> record(meta: meta + [tool: 'genomad'], table: table) },
            combined_metadata
        )

        vis_genomad = VIS_GENOMAD (
            tax_genomad
        )

        genomad_taxonomy = genomad_concat.csv.map { meta, table -> record(meta: meta, genomad: table) }
            .join(tax_genomad.map { r -> record(meta: withoutTool(r.meta), genomad_counts: r.taxonomy_counts) }, by: 'meta')
            .join(vis_genomad.map { r -> record(meta: withoutTool(r.meta), genomad_plots: r) }, by: 'meta')
            .map { r -> record(meta: r.meta, genomad: r.genomad, genomad_counts: r.genomad_counts, genomad_plots: r.genomad_plots) }
    }

    //
    // Taxonomy Viphogs
    //
    if ( !opts.skip_viphogs ) {
        viphogs = VIPHOGS_ANNOTATION (
            reps.map { r -> record(meta: r.meta, faa: r.faa, gff: r.reps_gff) },
            opts
        )

        tax_viphogs = TAX_VIPHOGS (
            viphogs.map { r -> record(meta: r.meta + [tool: 'viphogs'], table: r.taxonomy_tsv) },
            combined_metadata
        )

        vis_viphogs = VIS_VIPHOGS (
            tax_viphogs
        )

        viphogs_taxonomy = viphogs
            .join(tax_viphogs.map { r -> record(meta: withoutTool(r.meta), viphogs_counts: r.taxonomy_counts) }, by: 'meta')
            .join(vis_viphogs.map { r -> record(meta: withoutTool(r.meta), viphogs_plots: r) }, by: 'meta')
            .map { r ->
                record(
                    meta: r.meta,
                    viphogs: r.taxonomy_tsv,
                    viphogs_modified: r.modified_tsv,
                    viphogs_informative: r.informative_tsv,
                    viphogs_annotation: r.annotation_tsv,
                    viphogs_counts: r.viphogs_counts,
                    viphogs_plots: r.viphogs_plots
                )
            }
    }

    if ( !opts.skip_vitap ) {
        //
        // Taxonomy VITAP
        //
        vitap = VITAP (
            fna_chunks,
            opts.vitap_db
        )

        vitap_concat = CONCATENATE_VITAP (
            vitap.map { r -> tuple(r.meta, r.best_lineages) }.groupTuple(),
            1
        )

        tax_vitap = TAX_VITAP (
            vitap_concat.file_out.map { meta, table -> record(meta: meta + [tool: 'vitap'], table: table) },
            combined_metadata
        )

        vis_vitap = VIS_VITAP (
           tax_vitap
        )

        vitap_best = vitap_concat.file_out.map { meta, table -> record(meta: meta, vitap_best: table) }
            .join(tax_vitap.map { r -> record(meta: withoutTool(r.meta), vitap_counts: r.taxonomy_counts) }, by: 'meta')
            .join(vis_vitap.map { r -> record(meta: withoutTool(r.meta), vitap_plots: r) }, by: 'meta')
            .map { r -> record(meta: r.meta, vitap_best: r.vitap_best, vitap_counts: r.vitap_counts, vitap_plots: r.vitap_plots) }
    }

    // Steps above are optional: their fields are null when skipped
    taxonomy = reps
        .map { r -> record(meta: r.meta) }
        .join(vitap_best, by: 'meta', remainder: true)
        .join(genomad_taxonomy, by: 'meta', remainder: true)
        .join(viphogs_taxonomy, by: 'meta', remainder: true)
        .map { r ->
            record(
                meta: r.meta,
                vitap_best: r.vitap_best, vitap_counts: r.vitap_counts, vitap_plots: r.vitap_plots,
                genomad: r.genomad, genomad_counts: r.genomad_counts, genomad_plots: r.genomad_plots,
                viphogs: r.viphogs, viphogs_modified: r.viphogs_modified, viphogs_informative: r.viphogs_informative,
                viphogs_annotation: r.viphogs_annotation, viphogs_counts: r.viphogs_counts, viphogs_plots: r.viphogs_plots
            )
        }

    emit:
    taxonomy: Channel<TaxonomyResult> = taxonomy
}
