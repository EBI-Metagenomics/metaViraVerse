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
    tuple val(meta), path(fna)
    path(rna_inc_db)

    output:
    tuple val(meta), path("rna_candidates.tsv"), emit: rna_candidates
    tuple val("${task.process}"), val('infernal'), eval("cmsearch -h | sed -n 's/^# INFERNAL \\([0-9.]*\\).*/\\1/p'"), emit: versions_infernal, topic: versions

    script:
    """
    cmsearch --cpu ${task.cpus} --tblout rna_candidates.tsv ${rna_inc_db} ${fna}
    """
}
