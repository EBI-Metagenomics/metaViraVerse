process ANNOTATION {
    /*
     * Generate tabular file with ViPhOG annotation results for 
     * proteins predicted in viral contigs
    */
    tag "${meta.id}"

    label 'process_single'

    container 'quay.io/microbiome-informatics/virify-python3:1.1'

    input:
    tuple val(meta), path(tab), path(gff)

    output:
    tuple val(meta), path("*_annotation.tsv"), emit: annotations

    script:
    """
    viral_contigs_annotation.py -o . -t ${tab} -g ${gff}
    """
}
