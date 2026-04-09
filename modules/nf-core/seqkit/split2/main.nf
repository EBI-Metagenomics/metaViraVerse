process SEQKIT_SPLIT2 {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.8.1--h9ee0642_0' :
        'biocontainers/seqkit:2.8.1--h9ee0642_0' }"

    input:
    tuple val(meta), path(assembly)
    val(length)
    val(size)

    output:
    tuple val(meta), path("**/*.gz"), emit: assembly
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def has_length = length && length != []
    def has_size   = size && size != []

    // Validate mutually exclusive requirement
    if (has_length && has_size) {
        error("Cannot use both --by-length (${length}) and --by-size (${size}) at the same time.")
    }

    // Validate at least one is required
    if (!has_length && !has_size) {
        error("Must provide either 'length' or 'size' parameter. Both cannot be empty.")
    }

     // We are also tweaking the prefix to prevent names like <assembly_id>.part_001.gz to be used
     // in favour of <assembly_id>_part_001.gz which is more file name parsing friendly
     // which helps when concatenating chunked post-processed fasta files, such as the results of interposcan

    def chunk_by_length = has_length ? "--by-length ${length} --by-length-prefix ${meta.id}_" : ""
    def chunk_by_size   = has_size ? "--by-size ${size} --by-size-prefix ${meta.id}_" : ""
    """
    seqkit \\
        split2 \\
        $args ${chunk_by_length} ${chunk_by_size} \\
        --threads $task.cpus \\
        ${assembly} \\
        --out-dir ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(echo \$(seqkit 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}