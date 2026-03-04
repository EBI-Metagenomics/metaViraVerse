process METACEREBERUS {

    label 'process_high'
    tag "${meta.id}"
    container "quay.io/microbiome-informatics/metacerberus:1.4.0_1"

    input:
    tuple val(meta), path(faa)
    path metacerberus_db

    output:
    path "versions.yml",                         emit: versions

    // slurm execution? --hydraMPP-slurm $SLURM_JOB_NODELIST

    script:
    """
    metacerberus.py \\
       --protein ${faa} \\
       --dir-out ${meta.id}_metacerberus \\
       --replace \\
       --hmm ALL \\
       --db-path ${metacerberus_db} \\
       --cpus 16

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metacerberus: \$(metacerberus.py --version 2>&1 | sed 's/MetaCerberus: version: //g')
    END_VERSIONS
    """
}
