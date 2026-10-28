nextflow.enable.types = true

// One category's sequences, as emitted by SEPARATE_SEQUENCES
record Sequences {
    fna: Path
    gff: Path?
    faa: Path?
}

process CHOOSE_SEQUENCES {

    label 'process_low'
    tag "combined"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.75':
        'quay.io/biocontainers/biopython:1.75' }"

    input:
    record(
        meta: Map,
        map_file: Path
    )
    virus: Sequences
    prophage: Sequences
    plasmid: Sequences
    quality: Path?   // CheckV quality summary, null when --skip_checkv
    rna_gff: Path?   // barrnap GFF, null when --skip_rrna_detection

    // the rename map and the concatenated CheckV table are both named combined.tsv
    stage:
    stageAs map_file, 'mapping.tsv'

    output:
    record(
        meta: meta,
        metadata: file("${meta.id}_metadata.tsv"),
        filtered_metadata: file("${meta.id}_filtered.tsv"),
        excluded_metadata: file("${meta.id}_excluded.tsv"),
        filtered_fna: file("${meta.id}_filtered.fna"),
        filtered_gff: file("${meta.id}_filtered.gff"),
        filtered_faa: file("${meta.id}_filtered.faa"),
        virus_fna: file("${meta.id}_virus_filtered.fna"),
        virus_gff: file("${meta.id}_virus_filtered.gff"),
        virus_faa: file("${meta.id}_virus_filtered.faa"),
        prophage_fna: file("${meta.id}_prophage_filtered.fna"),
        prophage_gff: file("${meta.id}_prophage_filtered.gff"),
        prophage_faa: file("${meta.id}_prophage_filtered.faa"),
        plasmid_fna: file("${meta.id}_plasmid_filtered.fna"),
        plasmid_gff: file("${meta.id}_plasmid_filtered.gff"),
        plasmid_faa: file("${meta.id}_plasmid_filtered.faa"),
        stats: file("${meta.id}_stats.json")
    )

    topic:
    tuple(task.process, 'python', eval('python --version 2>&1 | sed "s/Python //g"')) >> 'versions'

    script:
    def virus_gff_arg    = virus.gff    ? "--viruses-gff ${virus.gff}"      : ""
    def virus_faa_arg    = virus.faa    ? "--viruses-faa ${virus.faa}"      : ""
    def prophage_gff_arg = prophage.gff ? "--prophages-gff ${prophage.gff}" : ""
    def prophage_faa_arg = prophage.faa ? "--prophages-faa ${prophage.faa}" : ""
    def plasmid_gff_arg  = plasmid.gff  ? "--plasmids-gff ${plasmid.gff}"   : ""
    def plasmid_faa_arg  = plasmid.faa  ? "--plasmids-faa ${plasmid.faa}"   : ""
    def rrna_arg    = rna_gff ? "--rrna ${rna_gff}"       : ""
    def quality_arg = quality ? "--quality ${quality}"    : ""
    def mapping_arg = "--map mapping.tsv"

    """
    choose_sequences.py \\
        --viruses ${virus.fna} \\
        ${virus_gff_arg} \\
        ${virus_faa_arg} \\
        --prophages ${prophage.fna} \\
        ${prophage_gff_arg} \\
        ${prophage_faa_arg} \\
        --plasmids ${plasmid.fna} \\
        ${plasmid_gff_arg} \\
        ${plasmid_faa_arg} \\
        ${rrna_arg} \\
        ${quality_arg} \\
        ${mapping_arg} \\
        --output-prefix ${meta.id}
    """
}
