process AMRINTEGRATOR {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:1.4.6--pyhdfd78af_0':
        'biocontainers/mgnify-pipelines-toolkit:1.4.6--pyhdfd78af_0' }"

    input:
    tuple val(meta),
     path(deeparg, stageAs: 'deeparg_hamr.tsv'),
     path(rgi, stageAs: 'rgi_hamr.tsv'),
     path(amrfp, stageAs: 'amrfinder.tsv'),
     path(gff)

    output:
    tuple val(meta), path("*.gff"), optional: true, emit: gff
    tuple val("${task.process}"), val('mgnify-pipelines-toolkit'), eval('get_mpt_version'), emit: versions_mgnify_pipelines_toolkit, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def deeparg_param = deeparg ? "--deeparg_hamr ${deeparg}" : ""
    def rgi_param = rgi ? "--rgi_hamr ${rgi}" : ""
    def amrfp_param = amrfp ? "--amrfp_out ${amrfp}" : ""
    """
    amr_integrator \\
        ${deeparg_param} \\
        ${rgi_param} \\
        ${amrfp_param} \\
        --cds_gff ${gff} \\
        --output integrated_${prefix}.gff \\
        $args
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args

    touch ${prefix}.gff
    """
}
