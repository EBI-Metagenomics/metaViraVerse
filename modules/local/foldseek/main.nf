/*
 * modules/local/foldseek/main.nf
 *
 * Structural homology search using Foldseek easy-search.
 * (van Kempen et al. Nature Methods 2024 doi:10.1038/s41592-023-02119-x)
 *
 * Two searches:
 *   1. BFVD (Big Fantastic Virus Database) — primary; always runs
 *      Kim et al. NAR 2024 doi:10.1093/nar/gkae1119
 *      Stream via bfvd.foldseek.com (no local download) or supply local DB path.
 *   2. PDB100 — optional; set params.foldseek_search_pdb = true
 *      Download once: foldseek databases PDB pdb_db tmp
 *
 * Output columns (both TSVs):
 *   query, target, evalue, bits, alntmscore, qtmscore, ttmscore,
 *   lddt, alnlen, pident
 *
 * TM-score interpretation: >0.5 = same fold; >0.7 = high confidence
 */

process FOLDSEEK_SEARCH {
    tag "${meta.id}"
    label 'process_high'

    conda 'bioconda::foldseek>=9.427df8a'
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/foldseek:9.427df8a--h4ac6f70_1' :
        'biocontainers/foldseek:9.427df8a--h4ac6f70_1' }"

    input:
    tuple val(meta), path(pdb_dir)
    path  bfvd_db

    output:
    tuple val(meta), path("${meta.id}_bfvd_hits.tsv"), emit: bfvd_hits
    tuple val(meta), path("${meta.id}_pdb_hits.tsv"),  emit: pdb_hits, optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    def threads  = task.cpus
    def run_pdb  = params.foldseek_search_pdb ?: false
    def pdb_db   = params.foldseek_pdb_db     ?: 'PDB'
    def evalue   = task.ext.evalue            ?: '0.001'
    """
    bash ${projectDir}/scripts/run_foldseek.sh \\
        "${pdb_dir}"               \\
        "${bfvd_db}"               \\
        "."                        \\
        "${threads}"               \\
        "${run_pdb}"               \\
        "${pdb_db}"                \\
        "${evalue}"

    mv bfvd_hits.tsv "${meta.id}_bfvd_hits.tsv"
    [ -f pdb_hits.tsv ] && mv pdb_hits.tsv "${meta.id}_pdb_hits.tsv" || true
    """

    stub:
    """
    printf "query\ttarget\tevalue\tbits\talntmscore\tqtmscore\tttmscore\tlddt\talnlen\tpident\n" \
        > "${meta.id}_bfvd_hits.tsv"
    printf "stub_seq.pdb\tBFVD_stub\t1e-10\t200\t0.82\t0.80\t0.84\t0.91\t120\t0.45\n" \
        >> "${meta.id}_bfvd_hits.tsv"
    touch "${meta.id}_pdb_hits.tsv"
    """
}
