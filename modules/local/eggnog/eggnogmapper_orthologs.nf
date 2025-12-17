process EGGNOGMAPPER_ORTHOLOGS {
    tag "${meta.id}"
    label 'process_long'

    conda "bioconda::eggnog-mapper=2.1.12"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/eggnog-mapper:2.1.12--pyhdfd78af_0':
        'biocontainers/eggnog-mapper:2.1.12--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)
    path(eggnog_data_dir)
    path(eggnog_db)
    path(eggnog_diamond_db)
    val(eggnog_database_version)

    output:
    tuple val(meta), path("*.hits")          , emit: hits
    tuple val(meta), path("*.seed_orthologs"), emit: orthologs
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def dbmem = task.memory.toMega() > 40000 ? '--dbmem' : ''

    def fasta_bool = fasta ? fasta.name : "no_fasta"
    def is_compressed = fasta_bool.endsWith(".gz")
    def fasta_name = fasta_bool.replace(".gz", "")
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi

    emapper.py \\
        ${args} \\
        -i ${fasta_name} \\
        --database ${eggnog_db} \\
        --dmnd_db ${eggnog_diamond_db} \\
        ${dbmem} \\
        --cpu ${task.cpus} \\
        --data_dir ${eggnog_data_dir} \\
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
    touch ${prefix}.emapper.seed_orthologs

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eggnog-mapper: \$(echo \$(emapper.py --version) | grep -o "emapper-[0-9]\\+\\.[0-9]\\+\\.[0-9]\\+" | sed "s/emapper-//")
        eggnog-mapper database: $eggnog_database_version
    END_VERSIONS
    """
}