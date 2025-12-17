process METACEREBERUS {

    label 'process_medium'
    tag "${meta.id}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/pip_metacerberus:304b52b16734d974':
        'community.wave.seqera.io/library/pip_metacerberus:5c70420d6a57973c' }"

    input:
    tuple val(meta), path(faa)

    output:
    tuple val(meta), path("*_renamed.fasta"),    emit: contigs_renamed
    path "versions.yml",                         emit: versions

    // slurm execution? --hydraMPP-slurm $SLURM_JOB_NODELIST

    script:
    """
    metacerberus.py \\
       --protein ${faa} \\
       --dir-out ${meta.id}_metacerberus \\
       --replace \\
       --hmm ALL \\
       --db-path ${params.metacerberus_db} \\
       --cpus 4

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metacerberus: \$(metacerberus.py --version 2>&1 | sed 's/MetaCerberus: version: //g')
    END_VERSIONS
    """
}
