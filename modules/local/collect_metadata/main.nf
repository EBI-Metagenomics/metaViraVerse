/*
 * module to generate viruses-all-metadata.tsv and plasmids-all-metadata.tsv
*/
process COLLECT_METADATA {

    label 'process_low'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.75':
        'quay.io/biocontainers/biopython:1.75' }"

    input:
    tuple val(meta_combined), path(combined_meta)
    tuple val(meta_viral), path(viral_clusters)
    tuple val(meta_plasmid), path(plasmid_clusters)
    path(additional_metadata)
    path(vitap_best)
    path(viphogs_taxonomy)
    path(genomad)
    path(map_file)
    path(iphop_genome)      // optional: [] when iPHoP was skipped
    path(iphop_genus)       // optional: [] when iPHoP was skipped
    path(spacepharer)       // optional: [] when SpacePHARER was not run

    output:
    tuple val(meta_viral), path("viruses-all-metadata.tsv.gz"),       emit: viruses_final_metadata
    tuple val(meta_plasmid), path("plasmids-all-metadata.tsv.gz"),    emit: plasmids_final_metadata
    tuple val(meta_plasmid), path("viruses-cluster-stats.tsv.gz"),    emit: viruses_reps_stats
    path "versions.yml",                       emit: versions

    script:
    def iphop_genome_arg = iphop_genome ? "--host-iphop-genome ${iphop_genome}" : ""
    def iphop_genus_arg  = iphop_genus  ? "--host-iphop-genus ${iphop_genus}"   : ""
    def spacepharer_arg  = spacepharer  ? "--host-spacepharer ${spacepharer}"   : ""
    def vitap_arg        = vitap_best  ? "--viruses_vitap ${vitap_best}"   : ""
    def viphogs_arg      = viphogs_taxonomy  ? "--viphogs_taxonomy ${viphogs_taxonomy}"   : ""
    def genomad_arg      = genomad  ? "--genomad ${genomad}"   : ""
    """
    collect_metadata.py \\
       --combined_meta ${combined_meta} \\
       --viruses_cluster ${viral_clusters} \\
       --plasmids_cluster ${plasmid_clusters} \\
       --additional_metadata ${additional_metadata} \\
       ${vitap_arg} \\
       ${viphogs_arg} \\
       ${genomad_arg} \\
       --map ${map_file} \\
       ${iphop_genome_arg} \\
       ${iphop_genus_arg} \\
       ${spacepharer_arg} \\
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
