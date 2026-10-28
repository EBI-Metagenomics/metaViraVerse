nextflow.enable.types = true

/*
 * module to generate viruses-all-metadata.tsv and plasmids-all-metadata.tsv
*/
process COLLECT_METADATA {

    label 'process_low'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.75':
        'quay.io/biocontainers/biopython:1.75' }"

    input:
    combined_meta: Path
    viral_meta: Map
    viral_clusters: Path
    plasmid_meta: Map
    plasmid_clusters: Path
    additional_metadata: Path
    vitap_best: Path?
    viphogs_taxonomy: Path?
    genomad: Path?
    map_file: Path
    iphop_genome: Path?      // null when iPHoP was skipped
    iphop_genus: Path?       // null when iPHoP was skipped
    spacepharer: Path?       // null when SpacePHARER was not run

    output:
    record(
        viral_meta: viral_meta,
        plasmid_meta: plasmid_meta,
        viruses_final_metadata: file("viruses-all-metadata.tsv.gz"),
        plasmids_final_metadata: file("plasmids-all-metadata.tsv.gz"),
        viruses_reps_stats: file("viruses-cluster-stats.tsv.gz")
    )

    topic:
    tuple(task.process, 'python', eval("python --version 2>&1 | sed 's/Python //g'")) >> 'versions'

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
    """
}
