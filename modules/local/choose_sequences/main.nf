process CHOOSE_SEQUENCES {

    label 'process_low'
    tag "combined"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.75':
        'quay.io/biocontainers/biopython:1.75' }"

    input:
    tuple val(meta_fna), path(fna_files)
    tuple val(meta_gff), path(gff_files)
    tuple val(meta_quality), path(quality, name: "quality_summary.tsv")
    tuple val(meta_rna_gff), path(rna_gff)
    tuple val(meta_map), path(map_file, name: "mapping.tsv")


    output:
    tuple val(meta_quality), path("${meta_quality.id}_metadata.tsv"),    emit: metadata
    tuple val(meta_quality), path("*filtered.fna"),                      emit: filtered_fna
    tuple val(meta_quality), path("*filtered.tsv"),                      emit: filtered_metadata
    tuple val(meta_quality), path("*filtered.gff"),                      emit: filtered_gff
    tuple val(meta_quality), path("*excluded.tsv"),                      emit: excluded_metadata
    path "versions.yml",                                                 emit: versions

    script:
    def fna_args   = fna_files.collect { it }.join(' ')
    def gff_args   = gff_files.collect { it }.join(' ')
    def rrna = rna_gff ? "--rrna ${rna_gff}" : ""
    def quality_arg = quality ? "--quality quality_summary.tsv" : ""
    def mapping = map_file ? "--map mapping.tsv" : ""

    """
    choose_sequences.py \\
        --fna ${fna_args} \\
        --gff ${gff_args} \\
        ${rrna} \\
        ${quality_arg} \\
        ${mapping} \\
        --output-prefix ${meta_quality.id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        biopython: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('biopython').version)")
    END_VERSIONS
    """
}
