nextflow.enable.types = true

include { CHECK_EXISTING_CATALOGUE } from '../../../modules/local/check_existing_catalogue'

/*
 * --mode update: read the catalogue_manifest.json of the existing catalogue
 * (--existing_catalogue, the --outdir of an earlier run), check its files and
 * hand them on as one record. The keys are the ones manifestInput() in
 * workflows/metaviraverse.nf writes.
*/

record ExistingCatalogue {
    dir: Path
    manifest: Path
    check_report: Path?           // existing_catalogue_check.json (set once CHECK_EXISTING_CATALOGUE passed)
    manifest_files: List<Path>    // every file the manifest lists (checked by CHECK_EXISTING_CATALOGUE)
    // accessions: an update continues after last_accession, with the same prefix and padding
    accession_prefix: String
    last_accession: Integer
    zero_padded_digits: Integer?
    // whole catalogue
    rename_map: Path
    catalogue_metadata: Path      // every unique sequence, with sequence_sha256
    excluded_metadata: Path?
    // viruses: vjdb-clustering tables (input of the clustering update)
    vjdb_clusters_csv: Path
    vjdb_reps_fna: Path
    // viruses: cluster representatives and their annotations
    reps_fna: Path
    reps_faa: Path
    reps_gff: Path
    annotated_gff: Path
    bacphlip: Path?
    genomad: Path?
    vitap: Path?
    viphogs: Path?
    viphogs_annotation: Path?
    iphop_genome: Path?
    iphop_genus: Path?
    spacepharer: Path?
    amrfinderplus: Path?
    deeparg: Path?
    rgi: Path?
    amr_gff: Path?
    phams: Path?
    hmmer_tables: List<Path>      // empty when --annotation_db was not set
    viruses_metadata: Path?
    viruses_cluster_stats: Path?
}

workflow LOAD_EXISTING_CATALOGUE {

    take:
    dir: Path
    rename_accession: String

    main:
    def loaded = loadExistingCatalogue(dir, rename_accession)
    log.info("--mode update: existing catalogue ${dir}, last accession ${loaded.accession_prefix}${loaded.last_accession}")

    check = CHECK_EXISTING_CATALOGUE(
        loaded.manifest,
        loaded.manifest_files,
        rename_accession
    )

    // downstream steps start only after the files were checked
    existing = check.map { r -> loaded + record(check_report: r.report) }

    emit:
    catalogue: Value<ExistingCatalogue> = existing
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Read the manifest, check what can be checked without reading the files (format,
// update readiness, accession prefix, required keys, files present) and build the record.
// Sizes and checksums are checked by CHECK_EXISTING_CATALOGUE.
def loadExistingCatalogue(dir: Path, rename_accession: String) -> ExistingCatalogue {
    def manifest_file = dir.resolve('catalogue_manifest.json')
    if( !manifest_file.exists() )
        error("--existing_catalogue ${dir} has no catalogue_manifest.json: give the --outdir of a metaViraVerse run that wrote one")

    def manifest = new groovy.json.JsonSlurper().parseText(manifest_file.text) as Map
    if( manifest.manifest_version != 1 )
        error("${manifest_file}: manifest_version ${manifest.manifest_version} is not supported (supported: 1)")

    def clustering = manifest.viral_clustering as Map
    if( !clustering.update_ready )
        error("${dir} cannot be updated: ${clustering.not_update_ready_reason}")

    def accessions = manifest.accessions as Map
    if( accessions.prefix != rename_accession )
        error("${dir} uses the accession prefix '${accessions.prefix}', but --rename_accession is '${rename_accession}': use the same prefix to update it")
    if( accessions.last == null )
        error("${manifest_file} has no accessions (empty catalogue)")

    def files = manifest.files as Map
    def required = [
        'rename_map', 'catalogue_metadata', 'viruses.vjdb_clusters_csv', 'viruses.vjdb_reps_fna',
        'viruses.reps_fna', 'viruses.reps_faa', 'viruses.reps_gff', 'viruses.annotated_gff'
    ]
    def missing_keys = required.findAll { k -> !files.containsKey(k) }
    if( missing_keys )
        error("${manifest_file} lists no ${missing_keys.join(', ')}")

    // a list key (viruses.hmmer_tables) holds several entries
    def all_files = files.values().collectMany { v -> (v instanceof List ? v : [v]) as List }
        .collect { e -> dir.resolve((e as Map).path as String) }
        .toList()
    def missing_files = all_files.findAll { f -> !f.exists() }
    if( missing_files )
        error("--existing_catalogue ${dir} is missing files listed in its manifest:\n  " + missing_files.join('\n  '))

    return record(
        dir: dir,
        manifest: manifest_file,
        check_report: null,
        manifest_files: all_files,
        accession_prefix: accessions.prefix as String,
        last_accession: accessions.last as Integer,
        zero_padded_digits: accessions.zero_padded_digits as Integer,
        rename_map: manifestPath(dir, files, 'rename_map'),
        catalogue_metadata: manifestPath(dir, files, 'catalogue_metadata'),
        excluded_metadata: manifestPath(dir, files, 'excluded_metadata'),
        vjdb_clusters_csv: manifestPath(dir, files, 'viruses.vjdb_clusters_csv'),
        vjdb_reps_fna: manifestPath(dir, files, 'viruses.vjdb_reps_fna'),
        reps_fna: manifestPath(dir, files, 'viruses.reps_fna'),
        reps_faa: manifestPath(dir, files, 'viruses.reps_faa'),
        reps_gff: manifestPath(dir, files, 'viruses.reps_gff'),
        annotated_gff: manifestPath(dir, files, 'viruses.annotated_gff'),
        bacphlip: manifestPath(dir, files, 'viruses.bacphlip'),
        genomad: manifestPath(dir, files, 'viruses.genomad'),
        vitap: manifestPath(dir, files, 'viruses.vitap'),
        viphogs: manifestPath(dir, files, 'viruses.viphogs'),
        viphogs_annotation: manifestPath(dir, files, 'viruses.viphogs_annotation'),
        iphop_genome: manifestPath(dir, files, 'viruses.iphop_genome'),
        iphop_genus: manifestPath(dir, files, 'viruses.iphop_genus'),
        spacepharer: manifestPath(dir, files, 'viruses.spacepharer'),
        amrfinderplus: manifestPath(dir, files, 'viruses.amrfinderplus'),
        deeparg: manifestPath(dir, files, 'viruses.deeparg'),
        rgi: manifestPath(dir, files, 'viruses.rgi'),
        amr_gff: manifestPath(dir, files, 'viruses.amr_gff'),
        phams: manifestPath(dir, files, 'viruses.phams'),
        hmmer_tables: ((files['viruses.hmmer_tables'] ?: []) as List)
            .collect { e -> dir.resolve((e as Map).path as String) }
            .toList(),
        viruses_metadata: manifestPath(dir, files, 'viruses.metadata'),
        viruses_cluster_stats: manifestPath(dir, files, 'viruses.cluster_stats')
    )
}

// Path of a single-file manifest key, or null when the step was skipped in that run
def manifestPath(dir: Path, files: Map, key: String) -> Path? {
    def entry = files[key] as Map
    return entry ? dir.resolve(entry.path as String) : null
}
