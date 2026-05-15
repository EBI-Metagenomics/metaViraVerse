/*
 * modules/local/ecod/main.nf
 *
 * ECOD domain annotation of predicted PDB structures.
 * Uses HHsearch against the ECOD HHM database fetched from Zenodo at runtime.
 *
 * ECOD Zenodo record (HHM DB): https://zenodo.org/records/13993145
 * Database is cached in workDir after first fetch — not re-downloaded between runs.
 *
 * Pre-filters structures to pLDDT > plddt_cutoff (default 0.7) before annotation.
 * SCOP cross-references are parsed from ECOD output internally.
 *
 * Cite: Cheng et al. NAR 2014
 */

process ECOD_ANNOTATE {
    tag "${meta.id}"
    label 'process_high'

    conda 'bioconda::hhsuite=3.3.0 conda-forge::python=3.10 conda-forge::requests'
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/hhsuite:3.3.0--py310pl5321h43eeafb_3' :
        'biocontainers/hhsuite:3.3.0--py310pl5321h43eeafb_3' }"

    input:
    tuple val(meta), path(pdb_dir)
    tuple val(meta), path(confidence_tsv)
    val   plddt_cutoff

    output:
    tuple val(meta), path("${meta.id}_ecod_annotations.tsv"), emit: ecod_tsv
    tuple val(meta), path("${meta.id}_scop_annotations.tsv"), emit: scop_tsv

    when:
    task.ext.when == null || task.ext.when

    script:
    def cutoff  = plddt_cutoff ?: 0.7
    def threads = task.cpus
    def zenodo  = "https://zenodo.org/records/13993145/files/ecod_hhm_db.tar.gz"
    """
    python3 ${projectDir}/scripts/annotate_ecod_hhsearch.py \\
        --pdb-dir "${pdb_dir}" \\
        --confidence-tsv "${confidence_tsv}" \\
        --plddt-cutoff ${cutoff} \\
        --output "${meta.id}_ecod_annotations.tsv" \\
        --scop-output "${meta.id}_scop_annotations.tsv" \\
        --zenodo-url "${zenodo}" \\
        --threads ${threads} \\
        --cache-dir "\${NXF_WORK:-\$PWD}/ecod_db_cache"
    """

    stub:
    """
    printf "seq_id\\tecod_uid\\tecod_xgroup\\tecod_hgroup\\tecod_tgroup\\tecod_fgroup\\tecod_evalue\\n" > "${meta.id}_ecod_annotations.tsv"
    printf "seq_id\\tscop_class\\tscop_fold\\tscop_superfamily\\tscop_family\\n" > "${meta.id}_scop_annotations.tsv"
    """
}
