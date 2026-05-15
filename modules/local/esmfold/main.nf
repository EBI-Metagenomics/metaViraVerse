/*
 * modules/local/esmfold/main.nf
 *
 * Protein structure prediction using ESMFold.
 * Two modes controlled by params.esmfold_mode:
 *   'api'   — ESM Metagenomic Atlas REST API (dev/small runs, max ~400aa)
 *   'local' — HuggingFace transformers local inference (production, up to ~1000aa)
 *
 * For local mode: label 'process_high_gpu' requires ~16GB VRAM (A100/V100).
 * ESMFold weights (~2.5 GB) downloaded once to HF_HOME at first run.
 */

process ESMFOLD {
    tag "${meta.id}"
    label params.esmfold_mode == 'local' ? 'process_high_gpu' : 'process_medium'

    conda 'conda-forge::python=3.10 conda-forge::requests'
    container "${ workflow.containerEngine == 'singularity' && params.esmfold_mode == 'local' ?
        'docker://ghcr.io/ebi-metagenomics/esmfold-local:latest' :
        'docker://ghcr.io/ebi-metagenomics/esmfold-api:latest' }"

    input:
    tuple val(meta), path(filtered_faa)

    output:
    tuple val(meta), path("${meta.id}_structures/"),      emit: pdb_dir
    tuple val(meta), path("${meta.id}_confidence.tsv"),   emit: confidence_tsv

    when:
    task.ext.when == null || task.ext.when

    script:
    def mode     = params.esmfold_mode ?: 'api'
    def delay    = task.ext.api_delay  ?: 1.2
    def script   = mode == 'local' ? 'run_esmfold_local.py' : 'run_esmfold_api.py'
    """
    mkdir -p "${meta.id}_structures"

    python3 ${projectDir}/scripts/${script} \\
        --input "${filtered_faa}" \\
        --outdir "${meta.id}_structures" \\
        --summary "${meta.id}_confidence.tsv" \\
        ${mode == 'api' ? "--delay ${delay}" : ""}
    """

    stub:
    """
    mkdir -p "${meta.id}_structures"
    touch "${meta.id}_structures/stub_seq.pdb"
    printf "seq_id\\tlength\\tmean_plddt\\tstatus\\n" > "${meta.id}_confidence.tsv"
    """
}
