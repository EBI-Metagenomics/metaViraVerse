process EGGNOGMAPPER_ANNOTATIONS {
    tag "${meta.id}"
    label 'process_long'

    conda "bioconda::eggnog-mapper=2.1.12"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/eggnog-mapper:2.1.12--pyhdfd78af_0':
        'biocontainers/eggnog-mapper:2.1.12--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(annotation_hit_table)
    path(eggnog_data_dir)
    val(eggnog_database_version)

    output:
    tuple val(meta), path("*.annotations")  , emit: annotations
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def dbmem = task.memory.toMega() > 40000 ? '--dbmem' : ''
    """
    # The --override option is used to prevent this process from failing on retries
    emapper.py \\
        ${args} \\
        --cpu ${task.cpus} ${dbmem} \\
        --data_dir ${eggnog_data_dir} \\
        --annotate_hits_table ${annotation_hit_table} \\
        --override \\
        --output ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eggnog-mapper: \$(echo \$(emapper.py --version) | grep -o "emapper-[0-9]\\+\\.[0-9]\\+\\.[0-9]\\+" | sed "s/emapper-//")
        eggnog-mapper database: $eggnog_database_version
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.emapper.annotations

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eggnog-mapper: \$(echo \$(emapper.py --version) | grep -o "emapper-[0-9]\\+\\.[0-9]\\+\\.[0-9]\\+" | sed "s/emapper-//")
        eggnog-mapper database: $eggnog_database_version
    END_VERSIONS
    """
}