process CHOOSE_SEQUENCES {

    label 'process_low'
    tag "combined"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.75':
        'quay.io/biocontainers/biopython:1.75' }"

    input:
    path(fna_files)
    path(gff_files)
    tuple val(meta_quality), path(quality)
    val(types)
    val(biomes)
    tuple val(meta_gff), path(rna_gff)
    path(map_file)


    output:
    tuple val(meta_quality), path("${meta_quality.id}.fna"),      emit: combined_fna
    tuple val(meta_quality), path("${meta_quality.id}.tsv"),      emit: metadata
    tuple val(meta_quality), path("filtered*.fna"),               emit: filtered_fna
    tuple val(meta_quality), path("filtered*.tsv"),               emit: filtered_metadata
    tuple val(meta_quality), path("${meta_quality.id}.gff"),      emit: filtered_gff
    path "versions.yml",                                          emit: versions

    script:
    def fna_args   = fna_files.collect { it }.join(' ')
    def gff_args   = gff_files.collect { it }.join(' ')
    def type_args  = types.join(' ')
    def biome_args = biomes.join(' ')
    def rrna = rna_gff ? "--rrna ${rna_gff}" : ""
    def quality_arg = quality ? "--quality ${quality}" : ""
    def mapping = map_file ? "--map ${map_file}" : ""

    """
    choose_sequences.py \\
        --fna ${fna_args} \\
        --gff ${gff_args} \\
        --type ${type_args} \\
        --biome ${biome_args} \\
        ${rrna} \\
        ${quality_arg} \\
        ${mapping} \\
        --output-fna ${meta_quality.id}.fna \\
        --output-tsv ${meta_quality.id}.tsv \\
        --output-gff ${meta_quality.id}.gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        biopython: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('biopython').version)")
    END_VERSIONS
    """
}
