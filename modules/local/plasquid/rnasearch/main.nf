nextflow.enable.types = true

process RNASEARCH {
    /*
     * Search plasmid contigs for Inc-associated RNA (Col/Rep) elements with Infernal
     * cmsearch, feeding both INCSEARCH's RNA-based Inc classification and (as a
     * sequence-context check) the protein-based Inc search.
    */

    label 'process_medium'
    tag "${meta.id}"
    container "docker://mgimenez720/plasquid:latest"

    input:
    record(
        meta: Map,
        fna: Path
    )
    rna_inc_db: Path

    output:
    record(
        meta: meta,
        rna_candidates: file("rna_candidates.tsv")
    )

    topic:
    tuple(task.process, 'infernal', eval("cmsearch -h | sed -n 's/^# INFERNAL \\([0-9.]*\\).*/\\1/p'")) >> 'versions'

    script:
    """
    cmsearch --cpu ${task.cpus} --tblout rna_candidates.tsv ${rna_inc_db} ${fna}
    """
}
