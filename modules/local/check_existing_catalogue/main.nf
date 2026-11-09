nextflow.enable.types = true

/*
 * --mode update: check that the files of the existing catalogue are the ones its
 * catalogue_manifest.json lists (size and SHA-256), before anything is built on them
*/
process CHECK_EXISTING_CATALOGUE {

    label 'process_single'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.75':
        'quay.io/biocontainers/biopython:1.75' }"

    input:
    manifest: Path
    files: List<Path>         // every file the manifest lists, under its published name
    rename_accession: String

    output:
    record(
        report: file("existing_catalogue_check.json")
    )

    script:
    """
    check_existing_catalogue.py \\
        --manifest ${manifest} \\
        --rename-accession ${rename_accession} \\
        -o existing_catalogue_check.json
    """

    stub:
    """
    echo '{"status": "ok", "manifest_sha256": "stub", "catalogue": {}}' > existing_catalogue_check.json
    """
}
