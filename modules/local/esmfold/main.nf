/*
 * modules/local/esmfold/main.nf
 *
 * Structure prediction using ESMFold (Lin et al. Science 2023).
 *
 * Two modes via params.esmfold_mode:
 *   'api'   — ESM Metagenomic Atlas REST API (dev/small runs; max ~400aa)
 *   'local' — HuggingFace transformers local inference (production; up to ~1000aa)
 *
 * Local mode requires ~16 GB VRAM; label 'process_high_gpu' is mandatory.
 * ESMFold weights (~2.5 GB) are downloaded once to HF_HOME at first run.
 * pLDDT values stored in B-factor column of ATOM records (0–100 scale).
 * Output confidence TSV normalises pLDDT to 0–1 for downstream filtering.
 *
 * Cite: Lin et al. Science 2023 doi:10.1126/science.ade2574
 */

process ESMFOLD {
    tag "${meta.id}"
    label params.esmfold_mode == 'local' ? 'process_high_gpu' : 'process_medium'

    conda (params.esmfold_mode == 'local'
        ? 'conda-forge::python=3.10 conda-forge::pytorch=2.0.1 huggingface::transformers=4.36.0 conda-forge::accelerate=0.24.0'
        : 'conda-forge::python=3.10 conda-forge::requests=2.31.0')

    container "${ workflow.containerEngine == 'singularity' ?
        ( params.esmfold_mode == 'local'
            ? 'docker://ghcr.io/ebi-metagenomics/esmfold-local:1.0.0'
            : 'docker://ghcr.io/ebi-metagenomics/esmfold-api:1.0.0' )
        : ( params.esmfold_mode == 'local'
            ? 'ghcr.io/ebi-metagenomics/esmfold-local:1.0.0'
            : 'ghcr.io/ebi-metagenomics/esmfold-api:1.0.0' ) }"

    input:
    tuple val(meta), path(filtered_faa)

    output:
    tuple val(meta), path("${meta.id}_structures/"),    emit: pdb_dir
    tuple val(meta), path("${meta.id}_confidence.tsv"), emit: confidence_tsv

    when:
    task.ext.when == null || task.ext.when

    script:
    def mode   = params.esmfold_mode ?: 'api'
    def delay  = task.ext.api_delay  ?: 1.2
    def script = (mode == 'local') ? 'run_esmfold_local.py' : 'run_esmfold_api.py'
    def extra  = (mode == 'api')   ? "--delay ${delay}"     : ''
    """
    mkdir -p "${meta.id}_structures"

    python3 ${projectDir}/scripts/${script} \\
        --input   "${filtered_faa}"         \\
        --outdir  "${meta.id}_structures"   \\
        --summary "${meta.id}_confidence.tsv" \\
        ${extra}
    """

    stub:
    """
    mkdir -p "${meta.id}_structures"
    # Minimal valid PDB with one ATOM record carrying a pLDDT B-factor of 85.0
    cat > "${meta.id}_structures/stub_seq.pdb" << 'PDBEOF'
ATOM      1  CA  MET A   1       1.000   1.000   1.000  1.00 85.00           C
END
PDBEOF
    printf "seq_id\tlength\tmean_plddt\tstatus\n"  > "${meta.id}_confidence.tsv"
    printf "stub_seq\t10\t0.850\tsuccess\n"        >> "${meta.id}_confidence.tsv"
    """
}
