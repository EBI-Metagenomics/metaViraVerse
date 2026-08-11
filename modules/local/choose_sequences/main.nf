process CHOOSE_SEQUENCES {

    label 'process_low'
    tag "combined"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.75':
        'quay.io/biocontainers/biopython:1.75' }"

    input:
    tuple val(meta_virus_fna),    path(virus_fna),    path(virus_gff),    path(virus_faa)
    tuple val(meta_prophage_fna), path(prophage_fna), path(prophage_gff), path(prophage_faa)
    tuple val(meta_plasmid_fna),  path(plasmid_fna),  path(plasmid_gff),  path(plasmid_faa)
    tuple val(meta_quality),      path(quality, name: "quality_summary.tsv")
    tuple val(meta_rna_gff),      path(rna_gff)
    tuple val(meta_map),          path(map_file, name: "mapping.tsv")

    output:
    tuple val(meta_map), path("${meta_map.id}_metadata.tsv"),          emit: metadata
    tuple val(meta_map), path("${meta_map.id}_filtered.tsv"),          emit: filtered_metadata
    tuple val(meta_map), path("${meta_map.id}_excluded.tsv"),          emit: excluded_metadata
    tuple val(meta_map), path("${meta_map.id}_filtered.fna"),          path("${meta_map.id}_filtered.gff"),          path("${meta_map.id}_filtered.faa"),          emit: filtered_data
    tuple val(meta_map), path("${meta_map.id}_virus_filtered.fna"),    path("${meta_map.id}_virus_filtered.gff"),    path("${meta_map.id}_virus_filtered.faa"),    emit: virus_data
    tuple val(meta_map), path("${meta_map.id}_prophage_filtered.fna"), path("${meta_map.id}_prophage_filtered.gff"), path("${meta_map.id}_prophage_filtered.faa"), emit: prophage_data
    tuple val(meta_map), path("${meta_map.id}_plasmid_filtered.fna"),  path("${meta_map.id}_plasmid_filtered.gff"),  path("${meta_map.id}_plasmid_filtered.faa"),  emit: plasmid_data
    tuple val("${task.process}"), val('python'), eval('python --version 2>&1 | sed "s/Python //g"'), topic: versions

    script:
    def virus_gff_arg    = virus_gff    ? "--viruses-gff ${virus_gff}"      : ""
    def virus_faa_arg    = virus_faa    ? "--viruses-faa ${virus_faa}"      : ""
    def prophage_gff_arg = prophage_gff ? "--prophages-gff ${prophage_gff}" : ""
    def prophage_faa_arg = prophage_faa ? "--prophages-faa ${prophage_faa}" : ""
    def plasmid_gff_arg  = plasmid_gff  ? "--plasmids-gff ${plasmid_gff}"   : ""
    def plasmid_faa_arg  = plasmid_faa  ? "--plasmids-faa ${plasmid_faa}"   : ""
    def rrna_arg    = rna_gff ? "--rrna ${rna_gff}"             : ""
    def quality_arg = quality ? "--quality quality_summary.tsv" : ""
    def mapping_arg = map_file ? "--map mapping.tsv"            : ""

    """
    choose_sequences.py \\
        --viruses ${virus_fna} \\
        ${virus_gff_arg} \\
        ${virus_faa_arg} \\
        --prophages ${prophage_fna} \\
        ${prophage_gff_arg} \\
        ${prophage_faa_arg} \\
        --plasmids ${plasmid_fna} \\
        ${plasmid_gff_arg} \\
        ${plasmid_faa_arg} \\
        ${rrna_arg} \\
        ${quality_arg} \\
        ${mapping_arg} \\
        --output-prefix ${meta_map.id}
    """
}
