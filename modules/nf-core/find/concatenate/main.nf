process FIND_CONCATENATE {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7f/7fd226561e12b32bcacdf4f5ff74577e76233adf52ae5cbc499a2cdfe0e27d82/data'
        : 'community.wave.seqera.io/library/findutils_pigz:c4dd5edc44402661'}"

    input:
    tuple val(meta), path(files_in, stageAs: 'to_concatenate/*', arity: '1..*')
    val number_of_header_lines

    output:
    tuple val(meta), path("${prefix}"), emit: file_out
    tuple val("${task.process}"), val("find"), eval("find --version | sed '1!d; s/.* //'"), topic: versions, emit: versions_find
    tuple val("${task.process}"), val("pigz"), eval("pigz --version 2>&1 | sed 's/pigz //g'"), topic: versions, emit: versions_pigz
    tuple val("${task.process}"), val("coreutils"), eval("cat --version | sed '1!d; s/.* //'"), topic: versions, emit: versions_coreutils

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""

    // | input     | output     | decomp      | compress at end |
    // |-----------|------------|-------------|-----------------|
    // | gzipped   | gzipped    | pigz -cd    | pigz            |
    // | ungzipped | ungzipped  | cat         |                 |
    // | gzipped   | ungzipped  | pigz -cd    |                 |
    // | ungzipped | gzipped    | cat         | pigz            |

    // Use input file ending as default
    // get file extensions, if extension is .gz then get the second to last extension as well
    file_extensions = files_in.collect { in_file -> in_file.name - in_file.getBaseName(in_file.name.endsWith('.gz') ? 2 : 1) }

    // Use input file ending as default for output file
    prefix = task.ext.prefix ?: "${meta.id}${file_extensions[0]}"

    if (files_in.any{ file -> file.toString().endsWith('.gz')} && !files_in.every{ file -> file.toString().endsWith('.gz') }) {
        error("All files provided to this module must either be gzipped (and have the .gz extension) or unzipped (and not have the .gz extension). A mix of both is not allowed.")
    }

    in_zip  = files_in[0].toString().endsWith('.gz')
    out_zip = task.ext.prefix ? task.ext.prefix.endsWith('.gz') : file_extensions[0].endsWith('.gz')

    // Always decompress for processing so header skipping works on plain text
    decomp   = in_zip ? "pigz -cd -p ${task.cpus}" : "cat"
    // Write to a plain (uncompressed) intermediate, then compress once at the end if needed
    out_bare = out_zip ? prefix.replaceAll(/\.gz$/, '') : prefix
    cmd_compress = out_zip ? "pigz -p ${task.cpus} ${args} ${out_bare}" : ""

    """
    first=true
    export header_lines=${number_of_header_lines}

    while IFS= read -r -d \$'\\0' file; do
        if [ "\$first" = true ]; then
            # First file: decompress and write everything (including header)
            ${decomp} \$file > ${out_bare}
            first=false
        else
            # Subsequent files: decompress and skip the first (header) line
            ${decomp} \$file | tail -n +\$((header_lines + 1)) >> ${out_bare}
        fi
    done < <( find to_concatenate/ -mindepth 1 -print0 | sort -z )

    ${cmd_compress}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"

    if (files_in.any{ file -> file.toString().endsWith('.gz')} && !files_in.every{ file -> file.toString().endsWith('.gz') }) {
        error("All files provided to this module must either be gzipped (and have the .gz extension) or unzipped (and not have the .gz extension). A mix of both is not allowed.")
    }

    """
    touch ${prefix}
    """
}
