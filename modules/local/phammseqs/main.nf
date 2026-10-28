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
    tuple val("${task.process}"), val('phammseqs'), val('1.3.0'), emit: versions_phammseqs, topic: versions
    tuple val("${task.process}"), val('clustalo'), eval("clustalo --help 2>&1 | grep 'Clustal Omega' | grep 'AndreaGiacomo' | sed 's/Clustal Omega - //g' | sed 's/(AndreaGiacomo)//g'"), emit: versions_clustalo, topic: versions
    tuple val("${task.process}"), val('mmseqs2'), eval("mmseqs -h 2>&1 | grep 'MMseqs2 Version' | sed 's/MMseqs2 Version: //g'"), emit: versions_mmseqs2, topic: versions

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
      -i ${meta.id}_pham/pham_fastas \\
      -o ${meta.id}_phams.tsv

    echo "Done."
    """
}
