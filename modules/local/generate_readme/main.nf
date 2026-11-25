nextflow.enable.types = true

/*
 * module to write README.md of the catalogue: summary numbers, tools that were run,
 * how the catalogue was built and links to the pipeline documentation
*/
process GENERATE_README {

    label 'process_single'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.75':
        'quay.io/biocontainers/biopython:1.75' }"

    input:
    record(
        catalogue_json: Path,
        viruses_metadata: Path,
        plasmids_metadata: Path,
        versions_yml: Path        // collated software versions: pipeline and tool versions
    )
    run_options: String           // JSON of the params that decide which tools ran
    catalogue_name: String?
    catalogue_version: String?
    catalogue_website_url: String?
    source_catalogues: String?    // comma-separated URLs

    output:
    record(
        readme: file("README.md")
    )

    // No `topic: ... >> 'versions'`: the input versions_yml is collated from that topic,
    // so emitting into it would make the pipeline wait forever (as for MultiQC).

    script:
    // single-quoted for the shell; a ' inside a value is written as '\''
    def name_arg    = catalogue_name        ? "--catalogue-name '${catalogue_name.replace("'", "'\\''")}'"               : ""
    def version_arg = catalogue_version     ? "--catalogue-version '${catalogue_version.replace("'", "'\\''")}'"         : ""
    def website_arg = catalogue_website_url ? "--website-url '${catalogue_website_url.replace("'", "'\\''")}'"           : ""
    def sources_arg = source_catalogues     ? "--source-catalogue '${source_catalogues.replace("'", "'\\''")}'"          : ""
    """
    cat > run_options.json <<'END_OPTIONS'
    ${run_options}
    END_OPTIONS

    generate_catalogue_readme.py \\
        --catalogue-json ${catalogue_json} \\
        --versions-yml ${versions_yml} \\
        --params run_options.json \\
        --metadata ${viruses_metadata} ${plasmids_metadata} \\
        ${name_arg} \\
        ${version_arg} \\
        ${website_arg} \\
        ${sources_arg} \\
        -o README.md
    """

    stub:
    """
    echo "# stub" > README.md
    """
}
