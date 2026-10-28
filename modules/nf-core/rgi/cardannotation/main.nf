process RGI_CARDANNOTATION {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/rgi:6.0.5--pyh05cac1d_0'
        : 'biocontainers/rgi:6.0.5--pyh05cac1d_0'}"

    input:
    path card

    output:
    path ("card_database_processed"), emit: db
    env 'RGI_VERSION', emit: tool_version
    env 'DB_VERSION', emit: db_version
    tuple val("${task.process}"), val('rgi'), env('RGI_VERSION'), emit: versions_rgi, topic: versions
    tuple val("${task.process}"), val('rgi-database'), env('DB_VERSION'), emit: versions_rgi_database, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    rgi card_annotation \\
        -i ${card}/card.json \\
        ${args}

    DB_VERSION=\$(ls card_database_*_all.fasta | sed "s/card_database_v\\([0-9].*[0-9]\\).*/\\1/")

    mkdir card_database_processed
    mv card*.fasta card_database_processed
    cp ${card}/* card_database_processed

    RGI_VERSION=\$(rgi main --version)
    """

    stub:
    """
    touch card.fasta
    touch card_all.fasta

    mkdir card_database_processed
    mv card*.fasta card_database_processed

    RGI_VERSION=\$(rgi main --version)
    DB_VERSION=stub_version
    """
}
