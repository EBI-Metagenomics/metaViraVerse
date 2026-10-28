nextflow.enable.types = true

// https://github.com/chg60/phammseqs
process PHAMMSEQS {

    label 'process_high'
    tag "${meta.id}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/clustalo_mmseqs2_pip_phammseqs:1b92ce4f391b1427':
        'community.wave.seqera.io/library/clustalo_mmseqs2_pip_phammseqs:147b2c0ddccc9182' }"

    input:
    record(
        meta: Map,
        faa: Path
    )

    output:
    record(
        meta: meta,
        phams_tsv: file("${meta.id}_phams.tsv"),
        phams_stats: file("${meta.id}_phammseqs.stdout")
    )

    topic:
    tuple(task.process, 'phammseqs', '1.3.0') >> 'versions'
    tuple(task.process, 'clustalo', eval("clustalo --help 2>&1 | grep 'Clustal Omega' | grep 'AndreaGiacomo' | sed 's/Clustal Omega - //g' | sed 's/(AndreaGiacomo)//g'")) >> 'versions'
    tuple(task.process, 'mmseqs2', eval("mmseqs -h 2>&1 | grep 'MMseqs2 Version' | sed 's/MMseqs2 Version: //g'")) >> 'versions'

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
