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
include { TAXONOMY_VISUALISATION as VIS_VIPHOGS        } from '../taxonomy_visualisation/main'
include { TAXONOMY_VISUALISATION as VIS_VITAP          } from '../taxonomy_visualisation/main'
include { TAXONOMY_VISUALISATION as VIS_GENOMAD        } from '../taxonomy_visualisation/main'



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

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
}

record GenomadTaxonomy {
    meta: Map
    genomad: Path
}

record ViphogsTaxonomy {
    meta: Map
    viphogs: Path
}

record TaxonomyResult {
    meta: Map
    vitap_best: Path?   // null when VITAP is skipped
    genomad: Path?      // null when geNomad is skipped
    viphogs: Path?      // null when ViPhOGs are skipped
}

workflow TAXONOMY_ASSIGNMENT {

    take:
    reps: Channel<TaxonomyInput>
    fna_chunks: Channel<FnaChunk>
    combined_metadata: Value<Path>
    skip_vitap: Boolean
    skip_genomad: Boolean
    skip_viphogs: Boolean
    viphog_db: Path?
    additional_viphog_model_data: Path?
    ncbi_db: Path?
    viphogs_factor_file: Path?

    main:

    def vitap_best: Channel<VitapBest> = channel.empty()
    def genomad_taxonomy: Channel<GenomadTaxonomy> = channel.empty()
    def viphogs_taxonomy: Channel<ViphogsTaxonomy> = channel.empty()

    def ch_fna_chunks = fna_chunks.map { r -> tuple(r.meta, r.fasta) }

    //
    // Taxonomy genoMAD
    //
    if ( !skip_genomad ) {
        genomad = GENOMAD_ENDTOEND (
            ch_fna_chunks,
            file(params.genomad_db)
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

        VIS_GENOMAD (
            tax_genomad
        )

        genomad_taxonomy = genomad_concat.csv.map { meta, table -> record(meta: meta, genomad: table) }
    }

    //
    // Taxonomy Viphogs
    //
    if ( !skip_viphogs ) {
        viphogs = VIPHOGS_ANNOTATION (
            reps.map { r -> record(meta: r.meta, faa: r.faa, gff: r.reps_gff) },
            viphog_db,
            additional_viphog_model_data,
            ncbi_db,
            viphogs_factor_file
        )

        tax_viphogs = TAX_VIPHOGS (
            viphogs.map { r -> record(meta: r.meta + [tool: 'viphogs'], table: r.taxonomy_tsv) },
            combined_metadata
        )

        VIS_VIPHOGS (
            tax_viphogs
        )

        viphogs_taxonomy = viphogs.map { r -> record(meta: r.meta, viphogs: r.taxonomy_tsv) }
    }

    if ( !skip_vitap ) {
        //
        // Taxonomy VITAP
        //
        vitap = VITAP (
            fna_chunks,
            file(params.vitap_db)
        )

        vitap_concat = CONCATENATE_VITAP (
            vitap.map { r -> tuple(r.meta, r.best_lineages) }.groupTuple(),
            1
        )

        tax_vitap = TAX_VITAP (
            vitap_concat.file_out.map { meta, table -> record(meta: meta + [tool: 'vitap'], table: table) },
            combined_metadata
        )

        VIS_VITAP (
           tax_vitap
        )

        vitap_best = vitap_concat.file_out.map { meta, table -> record(meta: meta, vitap_best: table) }
    }

    // Steps above are optional: their fields are null when skipped
    taxonomy = reps
        .map { r -> record(meta: r.meta) }
        .join(vitap_best, by: 'meta', remainder: true)
        .join(genomad_taxonomy, by: 'meta', remainder: true)
        .join(viphogs_taxonomy, by: 'meta', remainder: true)
        .map { r -> record(meta: r.meta, vitap_best: r.vitap_best, genomad: r.genomad, viphogs: r.viphogs) }

    emit:
    taxonomy: Channel<TaxonomyResult> = taxonomy
}
