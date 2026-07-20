/*
 * module to generate viruses-all-metadata.tsv and plasmids-all-metadata.tsv
*/
process COLLECT_METADATA {

    label 'process_low'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.75':
        'quay.io/biocontainers/biopython:1.75' }"

    input:
    tuple val(meta1), path(combined_meta)
    tuple val(meta2), path(viral_clusters)
    tuple val(meta3), path(plasmid_clusters)
    path(additional_metadata)
    tuple val(meta5), path(vitap_best)
    tuple val(meta_viphogs), path(viphogs_taxonomy)
    tuple val(meta_genomad), path(genomad)
    path(map_file)

    output:
    path("viruses-all-metadata.tsv.gz"),       emit: viruses_final_metadata
    path("plasmids-all-metadata.tsv.gz"),      emit: plasmids_final_metadata
    path("viruses-cluster-stats.tsv.gz"),      emit: viruses_reps_stats
    path "versions.yml",                       emit: versions

    script:
    """
    collect_metadata.py \\
       --combined_meta ${combined_meta} \\
       --viruses_cluster ${viral_clusters} \\
       --plasmids_cluster ${plasmid_clusters} \\
       --additional_metadata ${additional_metadata} \\
       --viruses_vitap ${vitap_best} \\
       --viphogs_taxonomy ${viphogs_taxonomy} \\
       --genomad ${genomad} \\
       --map ${map_file} \\
       --output_viruses viruses-all-metadata.tsv \\
       --output_plasmids plasmids-all-metadata.tsv \\
       --output_reps viruses-cluster-stats.tsv \\
       --compress

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
