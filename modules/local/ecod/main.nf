/*
 * modules/local/ecod/main.nf
 *
 * ECOD domain annotation of predicted PDB structures.
 * (Cheng et al. NAR 2014 doi:10.1093/nar/gkt1248)
 *
 * Uses HHsearch against the ECOD HHM database fetched from Zenodo at runtime.
 * Zenodo record: https://zenodo.org/records/13993145
 *   - Downloaded once; cached in params.ecod_db_cache (default: workDir/ecod_cache)
 *   - Not re-downloaded between runs if cache marker present
 *
 * Pre-filters structures to mean pLDDT >= params.plddt_cutoff (default 0.7)
 * before annotation. Structures below threshold written as 'low_plddt' — not
 * omitted — so the merge step sees every sequence.
 *
 * SCOP cross-references parsed from ECOD hit annotations at no extra cost.
 *
 * Cite: Cheng et al. Nucleic Acids Research 2014 doi:10.1093/nar/gkt1248
 */

process ECOD_ANNOTATE {
    tag "${meta.id}"
    label 'process_high'

    conda 'bioconda::hhsuite=3.3.0 conda-forge::python=3.10 conda-forge::requests=2.31.0'
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
    def cutoff    = plddt_cutoff                                        ?: 0.7
    def threads   = task.cpus
    def cache_dir = params.ecod_db_cache                                ?: "${workDir}/ecod_cache"
    def zenodo    = "https://zenodo.org/records/13993145/files/ecod_hhm_db.tar.gz"
    """
    python3 ${projectDir}/scripts/annotate_ecod_hhsearch.py \\
        --pdb-dir        "${pdb_dir}"                        \\
        --confidence-tsv "${confidence_tsv}"                 \\
        --plddt-cutoff   ${cutoff}                           \\
        --output         "${meta.id}_ecod_annotations.tsv"   \\
        --scop-output    "${meta.id}_scop_annotations.tsv"   \\
        --zenodo-url     "${zenodo}"                         \\
        --threads        ${threads}                          \\
        --cache-dir      "${cache_dir}"
    """

    stub:
    """
    printf "seq_id\tecod_uid\tecod_xgroup\tecod_hgroup\tecod_tgroup\tecod_fgroup\tecod_evalue\tstatus\n" \
        > "${meta.id}_ecod_annotations.tsv"
    printf "stub_seq\tECOD001\t2004.1\t2004.1.1\t2004.1.1.1\t2004.1.1.1.1\t1e-15\tannotated\n" \
        >> "${meta.id}_ecod_annotations.tsv"
    printf "seq_id\tscop_class\tscop_fold\tscop_superfamily\tscop_family\n" \
        > "${meta.id}_scop_annotations.tsv"
    printf "stub_seq\td\td.1\td.1.1\td.1.1.1\n" \
        >> "${meta.id}_scop_annotations.tsv"
    """
}
