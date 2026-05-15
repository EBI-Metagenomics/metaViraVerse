/*
 * modules/local/foldseek/main.nf
 *
 * Structural homology search using Foldseek easy-search.
 * Searches against:
 *   1. BFVD (Big Fantastic Virus Database) — primary viral structure DB
 *   2. PDB100 — experimental structures (optional toggle)
 *
 * BFVD streams via bfvd.foldseek.com at search time (no local download needed).
 * For offline/HPC use, pre-download and set params.bfvd_db to local path.
 *
 * Cite: Kim et al. NAR 2024 (doi:10.1093/nar/gkae1119)
 *       van Kempen et al. Nature Methods 2024
 */

process FOLDSEEK_SEARCH {
    tag "${meta.id}"
    label 'process_high'

    conda 'bioconda::foldseek>=9'
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/foldseek:9.427df8a--h4ac6f70_1' :
        'biocontainers/foldseek:9.427df8a--h4ac6f70_1' }"

    input:
    tuple val(meta), path(pdb_dir)
    path  bfvd_db   // path to BFVD foldseek DB directory; set to 'BFVD' to stream

    output:
    tuple val(meta), path("${meta.id}_bfvd_hits.tsv"),  emit: bfvd_hits
    tuple val(meta), path("${meta.id}_pdb_hits.tsv"),   emit: pdb_hits, optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    def threads     = task.cpus
    def run_pdb     = params.foldseek_search_pdb ?: false
    def pdb_db      = params.foldseek_pdb_db     ?: 'PDB'
    def fmt         = "query,target,evalue,bits,alntmscore,qtmscore,ttmscore,lddt,alnlen,pident"
    """
    bash ${projectDir}/scripts/run_foldseek.sh \\
        "${pdb_dir}" \\
        "${bfvd_db}" \\
        "." \\
        "${threads}" \\
        "${run_pdb}" \\
        "${pdb_db}"

    mv bfvd_hits.tsv "${meta.id}_bfvd_hits.tsv"
    [ -f pdb_hits.tsv ] && mv pdb_hits.tsv "${meta.id}_pdb_hits.tsv" || true
    """

    stub:
    """
    touch "${meta.id}_bfvd_hits.tsv"
    touch "${meta.id}_pdb_hits.tsv"
    """
}
