process RGI_DOWNLOADDB {
    label 'process_single'

    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://depot.galaxyproject.org/singularity/gnu-wget:1.18--h36e9172_9'
        : 'biocontainers/gnu-wget:1.18--h36e9172_9'}"

    output:
    path "card_dir/"   , emit: card_json
    tuple val("${task.process}"), val('wget'), eval("wget --version | head -1 | cut -d ' ' -f 3"), emit: versions_wget, topic: versions
    tuple val("${task.process}"), val('untar'), eval("echo \$(tar --version 2>&1) | sed 's/^.*(GNU tar) //; s/ Copyright.*\$//'"), emit: versions_untar, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir card_dir
    wget https://card.mcmaster.ca/latest/data
    tar -xvf data ./card.json
    mv card.json card_dir

    """
    stub:
    """
    mkdir card_dir
    touch card_dir/card.json
    """
}
