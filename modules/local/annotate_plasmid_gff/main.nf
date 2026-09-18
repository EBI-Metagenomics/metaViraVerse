process ANNOTATE_PLASMID_GFF {
    /*
     * Add plasmid-evidence attributes to a representative GFF:
     *  - RIP_domain/MOB_group/Inc_group (per-protein, from PLASQUID_WORKFLOW's
     *    protein_report.tsv) on matching CDS records
     *  - mob_suite_biomarker/mobsuite_identifier (per-contig, from MOB-suite's
     *    plasmids_biomarker_report.txt, optional) on matching sequence-level records
     *  - any attribute AMRINTEGRATOR added on top of this same GFF (AMR gene name,
     *    tool, coordinates, etc. -- computed as "keys present in --amr-gff but not
     *    in --gff", optional) on matching records of any feature type
    */

    label 'process_low'
    tag "${meta.id}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.75':
        'quay.io/biocontainers/biopython:1.75' }"

    input:
    tuple val(meta), path(gff), path(protein_report), path(biomarker_report), path(amr_gff)

    output:
    tuple val(meta), path("${meta.id}_plasquid_annotated.gff"), emit: gff
    path "versions.yml",                                        emit: versions

    script:
    def biomarker_report_arg = biomarker_report ? "--biomarker-report ${biomarker_report}" : ""
    def amr_gff_arg          = amr_gff          ? "--amr-gff ${amr_gff}"                   : ""
    """
    annotate_plasmid_gff.py \\
        --gff ${gff} \\
        --table ${protein_report} \\
        ${biomarker_report_arg} \\
        ${amr_gff_arg} \\
        --output ${meta.id}_plasquid_annotated.gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
