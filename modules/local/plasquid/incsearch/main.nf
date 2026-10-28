nextflow.enable.types = true

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
    record(
        meta: Map,
        proteins: Path,
        rna_candidates: Path,
        protein_to_contig: Path
    )
    incsearch_db: Path

    output:
    record(
        meta: meta,
        inc_candidates: file("inc_candidates.tsv"),
        inc_classification: file("classification_table.tsv"),
        filt_classification: file("filtered_classif.tsv")
    )

    topic:
    tuple(task.process, 'hmmer', eval("hmmsearch -h | sed -n 's/^# HMMER \\([0-9.]*\\).*/\\1/p'")) >> 'versions'
    tuple(task.process, 'r-base', eval("R --version | sed -n '1s/^R version \\([0-9.]*\\).*/\\1/p'")) >> 'versions'

    script:
    """
    hmmsearch -o log --cpu ${task.cpus} --domtblout inc_candidates.tsv ${incsearch_db} ${proteins}

    plasquid_inc_classification.R inc_candidates.tsv ${rna_candidates} ${protein_to_contig}
    plasquid_filter_classification.R classification_table.tsv
    """
}
