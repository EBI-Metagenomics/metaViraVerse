process CRISPRCAS_FINDER {

    tag "${meta.id}"

    container 'quay.io/microbiome-informatics/genomes-pipeline.crisprcasfinder:4.3.2'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("crisprcasfinder_results/${meta.id}_crisprcasfinder.gff"),    emit: gff
    tuple val(meta), path("crisprcasfinder_results/${meta.id}_crisprcasfinder.tsv"),    emit: tsv
    tuple val(meta), path("crisprcasfinder_results/${meta.id}_crisprcasfinder_hq.gff"), emit: hq_gff
    path "versions.yml",                                                                       emit: versions

    script:
    """
    # Remove results folder if it already exists to prevent restarts from failing
    if [ -d "crisprcasfinder_results" ]; then
        rm -rf "crisprcasfinder_results"
    fi

    CRISPRCasFinder.pl -i ${fasta} \
    -so /opt/CRISPRCasFinder/sel392v2.so \
    -def G \
    -drpt /opt/CRISPRCasFinder/supplementary_files/repeatDirection.tsv \
    -outdir crisprcasfinder_results

    echo "Running post-processing"

    process_crispr_results.py \
    --tsv-report crisprcasfinder_results/TSV/Crisprs_REPORT.tsv \
    --gffs crisprcasfinder_results/GFF/*gff \
    --tsv-output crisprcasfinder_results/${meta.id}_crisprcasfinder.tsv \
    --gff-output crisprcasfinder_results/${meta.id}_crisprcasfinder.gff \
    --gff-output-hq crisprcasfinder_results/${meta.id}_crisprcasfinder_hq.gff \
    --fasta $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CRISPRCasFinder: \$(CRISPRCasFinder.pl -v 2>&1 | grep -oE '[0-9]+\\.[0-9]+\\.[0-9]+')
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS

    """
}
