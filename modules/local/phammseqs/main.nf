// https://github.com/chg60/phammseqs
process PHAMMSEQS {

    label 'process_high'
    tag "${meta.id}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/clustalo_mmseqs2_pip_phammseqs:1b92ce4f391b1427':
        'community.wave.seqera.io/library/clustalo_mmseqs2_pip_phammseqs:147b2c0ddccc9182' }"

    input:
    tuple val(meta), path(faa)

    output:
    tuple val(meta), path("${meta.id}_phams.tsv"),        emit: phams_tsv
    tuple val(meta), path("${meta.id}_phammseqs.stdout"), emit: phams_stats
    path "versions.yml",                                  emit: versions

    script:
    """
    # Run with -p (required to modify input) for file strain_genes.tsv or gene_presence_absence.csv
    echo "Running phammseqs"
    phammseqs -v \\
      --outdir ${meta.id}_pham \\
      --cluster-mode 2 \\
      ${faa} > ${meta.id}_phammseqs.stdout

    echo "Running phams_to_tsv"
    phams_to_tsv.py \\
      -i ${meta.id}_pham \\
      -o ${meta.id}_phams.tsv

    echo "Done."

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phammseqs: 1.3.0
        clustalo: \$(clustalo --help 2>&1 | grep 'Clustal Omega' | grep 'AndreaGiacomo' | sed 's/Clustal Omega - //g' | sed 's/(AndreaGiacomo)//g')
        mmseqs2: \$(mmseqs -h 2>&1 | grep 'MMseqs2 Version' | sed 's/MMseqs2 Version: //g')
    END_VERSIONS
    """
}
