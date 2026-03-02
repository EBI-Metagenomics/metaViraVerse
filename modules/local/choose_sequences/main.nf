process CHOOSE_SEQUENCES {

    label 'process_low'
    tag "combined"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.75':
        'quay.io/biocontainers/biopython:1.75' }"

    input:
    path(fna_files)
    val(types)
    val(biomes)

    output:
    path("combined.fna"),      emit: combined_fna
    path("combined_meta.tsv"), emit: metadata
    path "versions.yml",       emit: versions

    script:
    def fna_args   = fna_files.collect { it }.join(' ')
    def type_args  = types.join(' ')
    def biome_args = biomes.join(' ')
    """
    choose_sequences.py \\
        --fna ${fna_args} \\
        --type ${type_args} \\
        --biome ${biome_args} \\
        --output-fna combined.fna \\
        --output-tsv combined_meta.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        biopython: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('biopython').version)")
    END_VERSIONS
    """
}
