nextflow.enable.types = true

/*
 * Write catalogue_manifest.json: the accession range, deduplication and clustering
 * settings, software versions, and the publish path, size and SHA-256 of every file
 * a later `--mode update` run reads from this catalogue.
 *
 * No `topic: versions` here: this process reads the collated versions, which are
 * complete only once every process that sends to the topic has finished.
*/
process WRITE_CATALOGUE_MANIFEST {

    label 'process_single'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.75':
        'quay.io/biocontainers/biopython:1.75' }"

    input:
    record(
        files: List<Path>,      // the published files (unique file names)
        entries: List<String>   // one "KIND KEY TARGET NAME" per file, see bin/write_catalogue_manifest.py
    )
    rename_accession: String
    viral_clustering: String
    mode: String              // --mode: 'generate' or 'update'
    parent_check: Path?       // --mode update: existing_catalogue_check.json of the catalogue this run updates

    output:
    record(
        catalogue_manifest: file("catalogue_manifest.json")
    )

    script:
    def entry_args = entries.collect { e -> "--entry ${e}" }.join(' \\\n        ')
    def parent_arg = parent_check ? "--parent-check ${parent_check}" : ""
    """
    write_catalogue_manifest.py \\
        --rename-accession ${rename_accession} \\
        --viral-clustering ${viral_clustering} \\
        --mode ${mode} \\
        ${parent_arg} \\
        --pipeline-version '${workflow.manifest.version}' \\
        --run-name '${workflow.runName}' \\
        --commit '${workflow.commitId ?: ''}' \\
        ${entry_args} \\
        -o catalogue_manifest.json
    """

    stub:
    """
    echo '{}' > catalogue_manifest.json
    """
}
