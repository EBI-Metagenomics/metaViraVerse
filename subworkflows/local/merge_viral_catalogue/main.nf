nextflow.enable.types = true

include { MERGE_UPDATE_DROP_IDS ; MERGE_UPDATE ; MERGE_UPDATE_REPS } from '../../../modules/local/merge_catalogue_update'

include { INDEX_RESULTS                    } from '../index_results/main'
include { ExistingCatalogue                } from '../existing_catalogue/main'
include { ViralCatalogue                   } from '../process_viral_sequences/main'

/*
 * --mode update: the viral catalogue of the update covers the whole catalogue in its
 * clustering, but its annotations only the representatives of new clusters. Merge each
 * per-representative file with the existing catalogue's (existing rows without the
 * representatives dropped by cluster merges, then the new rows) and index the merged
 * sequences and GFF again.
 *
 * Outputs computed from the whole set (taxonomy counts, Krona/Sankey/iTOL plots, host
 * Sankey plots, HMMER summaries, phams, ViPhOGs intermediate tables) are not merged
 * yet: they are left out (null) instead of being published for the new clusters only.
*/

record MergeJob {
    meta: Map
    key: String
    format: String
    id_column: String?
    id_index: Integer?
    existing: Path
    update: Path
    name: String
}

record MergedFile {
    meta: Map
    key: String
    merged: Path
}

workflow MERGE_VIRAL_CATALOGUE {

    take:
    catalogue: Channel<ViralCatalogue>    // of the update: annotations of the new clusters only
    existing: Value<ExistingCatalogue?>

    main:

    // one viral catalogue: a dataflow value, so that every merge job reads the same inputs
    def update = catalogue.combine(existing: existing).reduce { a, _b -> a }

    // representatives absorbed by merges and their proteins: their rows are dropped
    def drop_ids = MERGE_UPDATE_DROP_IDS(
        update.map { r -> r.meta },
        update.map { r -> r.dropped_reps },
        update.map { r -> r.existing.reps_gff }
    )

    def merged = MERGE_UPDATE(
        update.flatMap { r -> mergeJobs(r as ViralCatalogue, r.existing) },
        drop_ids
    )

    def reps_list = MERGE_UPDATE_REPS(
        update.map { r -> r.meta },
        update.map { r -> r.clusters_tsv }
    )

    // a channel again, to combine with the merged files
    def files = update
        .flatMap { r -> [r] }
        .combine(merged_files: merged.collect(), reps_list: reps_list)

    // index the merged representatives as in generate mode
    def indexed = INDEX_RESULTS(
        files.map { r ->
            record(
                meta: r.meta,
                reps_gff: mergedFile(r.merged_files, 'reps_gff'),
                fna: mergedFile(r.merged_files, 'fna'),
                faa: mergedFile(r.merged_files, 'faa')
            )
        },
        files.map { r -> record(meta: r.meta, bacphlip: mergedFile(r.merged_files, 'bacphlip')) }
    )

    def merged_catalogue = files
        .join(indexed.map { r -> record(meta: r.meta, index: r) }, by: 'meta')
        .map { r ->
            def m = r.merged_files
            record(
                meta: r.meta,
                clusters_tsv: r.clusters_tsv, annotate_tsv: r.annotate_tsv,
                ani_tsv: r.ani_tsv, ani_clusters_tsv: r.ani_clusters_tsv,
                vjdb_clusters_csv: r.vjdb_clusters_csv, vjdb_reps_fna: r.vjdb_reps_fna,
                vjdb_new_clusters: r.vjdb_new_clusters, cluster_events: r.cluster_events,
                dropped_reps: r.dropped_reps,
                reps_list: r.reps_list,
                fna_gz: r.index.fna_gz,
                faa_gz: r.index.faa_gz,
                reps_gff: mergedFile(m, 'reps_gff'),
                final_gff: mergedFile(m, 'annotated_gff'),
                index: r.index,
                taxonomy: record(
                    meta: r.meta,
                    vitap_best: mergedFile(m, 'vitap'), vitap_counts: null, vitap_plots: null,
                    genomad: mergedFile(m, 'genomad'), genomad_counts: null, genomad_plots: null,
                    viphogs: mergedFile(m, 'viphogs'), viphogs_modified: null, viphogs_informative: null,
                    viphogs_annotation: mergedFile(m, 'viphogs_annotation'), viphogs_counts: null, viphogs_plots: null
                ),
                hosts: record(
                    meta: r.meta,
                    iphop_genome: mergedFile(m, 'iphop_genome'), iphop_genome_sankey: null, iphop_genome_reduced_sankey: null,
                    iphop_genus: mergedFile(m, 'iphop_genus'), iphop_genus_sankey: null, iphop_genus_reduced_sankey: null,
                    spacepharer: mergedFile(m, 'spacepharer'), spacepharer_predictions: null, spacepharer_lineage_json: null,
                    spacepharer_sankey: null, spacepharer_reduced_sankey: null
                ),
                proteins: record(
                    meta: r.meta,
                    hmmer_tables: r.proteins.hmmer_tables != null ? mergedFiles(m, 'hmmer:') : null,
                    hmmer_summaries: null,
                    amr_gff: mergedFile(m, 'amr_gff'),
                    amrfinderplus_tsv: mergedFile(m, 'amrfinderplus'),
                    deeparg_tsv: mergedFile(m, 'deeparg'),
                    rgi_tsv: mergedFile(m, 'rgi'),
                    phams_tsv: null, phams_plots: null,
                    mmseqs_tsv: null, mmseqs_plots: null
                )
            )
        }

    emit:
    catalogue: Channel<ViralCatalogue> = merged_catalogue
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Name of a merged file: the update's file name, uncompressed for the files INDEX_RESULTS compresses
def plainName(f: Path) -> String {
    def name = f.name as String
    return name.endsWith('.gz') ? name.substring(0, name.length() - 3) : name
}

// One merge job, or none when neither catalogue has the file. A file in only one of
// them means that the update was annotated with other tools than the existing catalogue.
def mergeJob(meta: Map, key: String, format: String, id_column: String?, id_index: Integer?, existing: Path?, update: Path?, name: String?) -> List<MergeJob> {
    if( existing == null && update == null )
        return []
    if( existing == null || update == null )
        error("--mode update: '${key}' is ${existing == null ? 'not in the existing catalogue' : 'missing from the update'}; " +
              "annotate the update with the same tools and databases as the existing catalogue")
    return [record(meta: meta, key: key, format: format, id_column: id_column, id_index: id_index,
                   existing: existing, update: update, name: name ?: (update.name as String))]
}

// The files of the update to merge with the existing catalogue, with the ID of their rows
// (see "File formats" in CLAUDE.md)
def mergeJobs(c: ViralCatalogue, e: ExistingCatalogue) -> List<MergeJob> {
    def m = c.meta
    def jobs = [
        mergeJob(m, 'reps_gff', 'gff', null, null, e.reps_gff, c.reps_gff, plainName(c.reps_gff)),
        mergeJob(m, 'fna', 'fasta', null, null, e.reps_fna, c.fna_gz, plainName(c.fna_gz)),
        mergeJob(m, 'faa', 'fasta', null, null, e.reps_faa, c.faa_gz, plainName(c.faa_gz)),
        mergeJob(m, 'annotated_gff', 'gff', null, null, e.annotated_gff, c.final_gff, null),
        mergeJob(m, 'bacphlip', 'tsv', null, 0, e.bacphlip, c.index.bacphlip_gz, c.index.bacphlip_gz ? plainName(c.index.bacphlip_gz) : null),
        mergeJob(m, 'genomad', 'tsv', 'seq_name', null, e.genomad, c.taxonomy.genomad, null),
        mergeJob(m, 'vitap', 'tsv', 'Genome_ID', null, e.vitap, c.taxonomy.vitap_best, null),
        mergeJob(m, 'viphogs', 'tsv', 'contig_ID', null, e.viphogs, c.taxonomy.viphogs, null),
        mergeJob(m, 'viphogs_annotation', 'tsv', 'Contig', null, e.viphogs_annotation, c.taxonomy.viphogs_annotation, null),
        mergeJob(m, 'iphop_genome', 'csv', 'Virus', null, e.iphop_genome, c.hosts.iphop_genome, null),
        mergeJob(m, 'iphop_genus', 'csv', 'Virus', null, e.iphop_genus, c.hosts.iphop_genus, null),
        mergeJob(m, 'spacepharer', 'tsv', 'viral_seq', null, e.spacepharer, c.hosts.spacepharer, null),
        mergeJob(m, 'amrfinderplus', 'tsv', 'Protein id', null, e.amrfinderplus, c.proteins.amrfinderplus_tsv, null),
        mergeJob(m, 'deeparg', 'tsv', 'input_sequence_id', null, e.deeparg, c.proteins.deeparg_tsv, null),
        mergeJob(m, 'rgi', 'tsv', 'input_sequence_id', null, e.rgi, c.proteins.rgi_tsv, null),
        mergeJob(m, 'amr_gff', 'gff', null, null, e.amr_gff, c.proteins.amr_gff, null)
    ]
    // HMMER: one table per HMM database, matched by file name
    def old_tables: List<Path> = e.hmmer_tables
    def new_tables: List<Path> = c.proteins.hmmer_tables ?: []
    def names = (old_tables + new_tables).collect { f -> f.name as String }.toSet()
    def hmmer = names.collect { n ->
        mergeJob(m, "hmmer:${n}", 'tbl', null, null, old_tables.find { f -> f.name == n }, new_tables.find { f -> f.name == n }, n)
    }
    return (jobs + hmmer.toList()).collectMany { j -> j }.toList()
}

def mergedFile(files: Bag<MergedFile>, key: String) -> Path? {
    return files.find { f -> f.key == key }?.merged
}

def mergedFiles(files: Bag<MergedFile>, prefix: String) -> List<Path> {
    return files.findAll { f -> f.key.startsWith(prefix) }.collect { f -> f.merged }.toList()
}
