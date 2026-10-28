process INCSEARCH {
    /*
     * Classify plasmid incompatibility (Inc) groups from two lines of evidence:
     * hmmsearch of predicted proteins against Inc-associated protein profiles, and
     * RNASEARCH's cmsearch hits against Inc-associated RNA covariance models. Each
     * Inc family has its own curated score cutoff, and known cross-reactive pairs
     * (IncFII/IncZ, Col/Inc13) are resolved by keeping only the primary call.
    */

    label 'process_medium'
    tag "${meta.id}"
    container "docker://mgimenez720/plasquid:latest"

    input:
    tuple val(meta), path(proteins), path(rna_candidates), path(protein_to_contig)
    path(incsearch_db)

    output:
    tuple val(meta), path("inc_candidates.tsv"),       emit: inc_candidates
    tuple val(meta), path("classification_table.tsv"), emit: inc_classification
    tuple val(meta), path("filtered_classif.tsv"),     emit: filt_classification
    tuple val("${task.process}"), val('hmmer'), eval("hmmsearch -h | sed -n 's/^# HMMER \\([0-9.]*\\).*/\\1/p'"), emit: versions_hmmer, topic: versions
    tuple val("${task.process}"), val('r-base'), eval("R --version | sed -n '1s/^R version \\([0-9.]*\\).*/\\1/p'"), emit: versions_r_base, topic: versions

    script:
    """
    hmmsearch -o log --cpu ${task.cpus} --domtblout inc_candidates.tsv ${incsearch_db} ${proteins}

    plasquid_inc_classification.R inc_candidates.tsv ${rna_candidates} ${protein_to_contig}
    plasquid_filter_classification.R classification_table.tsv
    """
}
