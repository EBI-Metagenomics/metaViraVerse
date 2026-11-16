nextflow.enable.types = true

/*
 * Test wrapper: nf-test starts a workflow with `run(inputs)` from an untyped script, which
 * leaves its inputs as plain dataflow objects, so the named `combine` of MERGE_VIRAL_CATALOGUE
 * fails on them. This typed wrapper calls LOAD_EXISTING_CATALOGUE and MERGE_VIRAL_CATALOGUE
 * as the pipeline does.
 */

include { LOAD_EXISTING_CATALOGUE } from '../../existing_catalogue/main'
include { MERGE_VIRAL_CATALOGUE   } from '../main'
include { ViralCatalogue          } from '../../process_viral_sequences/main'

workflow TEST_MERGE_VIRAL_CATALOGUE {

    take:
    existing_dir: Path   // the existing catalogue (--existing_catalogue)
    data: Path           // the update's viral catalogue (see tests/data/README)

    main:
    def meta = [id: 'viruses']
    def update = channel.of(record(
        meta: meta,
        clusters_tsv: data.resolve('viruses_clusters.tsv'),
        annotate_tsv: data.resolve('viruses_clusters.tsv'),
        ani_tsv: null, ani_clusters_tsv: null,
        vjdb_clusters_csv: null, vjdb_reps_fna: null, vjdb_new_clusters: null, cluster_events: null,
        dropped_reps: data.resolve('viruses_dropped_reps.txt'),
        reps_list: data.resolve('viruses_clusters.tsv'),
        fna_gz: data.resolve('_viruses_reps.fa.gz'),
        faa_gz: data.resolve('_viruses_reps.faa.gz'),
        reps_gff: data.resolve('viruses_reps.gff'),
        final_gff: data.resolve('viruses_reps_annotated.gff.gz'),
        // only the fields MERGE_VIRAL_CATALOGUE reads are set; the others are null or a placeholder
        index: record(meta: meta,
            fna_gz: data.resolve('_viruses_reps.fa.gz'), fna_fai: null, fna_gzi: null,
            faa_gz: data.resolve('_viruses_reps.faa.gz'), faa_fai: null, faa_gzi: null,
            gff_gz: data.resolve('viruses_reps_annotated.gff.gz'), gff_tbi: data.resolve('README'),
            bacphlip_gz: data.resolve('viruses.bacphlip.gz'), bacphlip_tbi: null),
        taxonomy: record(meta: meta,
            vitap_best: null, vitap_counts: null, vitap_plots: null,
            genomad: null, genomad_counts: null, genomad_plots: null,
            viphogs: null, viphogs_modified: null, viphogs_informative: null,
            viphogs_annotation: null, viphogs_counts: null, viphogs_plots: null),
        hosts: record(meta: meta,
            iphop_genome: null, iphop_genome_sankey: null, iphop_genome_reduced_sankey: null,
            iphop_genus: null, iphop_genus_sankey: null, iphop_genus_reduced_sankey: null,
            spacepharer: null, spacepharer_predictions: null, spacepharer_lineage_json: null,
            spacepharer_sankey: null, spacepharer_reduced_sankey: null),
        proteins: record(meta: meta,
            hmmer_tables: null, hmmer_summaries: null, amr_gff: null,
            amrfinderplus_tsv: null, deeparg_tsv: null, rgi_tsv: null,
            phams_tsv: null, phams_plots: null, mmseqs_tsv: null, mmseqs_plots: null)
    ))
    // annotations other than BACPHLIP were skipped in both catalogues
    def existing = LOAD_EXISTING_CATALOGUE(existing_dir, 'MGYV', [:])
    def merged = MERGE_VIRAL_CATALOGUE(update, existing)

    emit:
    catalogue: Channel<ViralCatalogue> = merged
}
