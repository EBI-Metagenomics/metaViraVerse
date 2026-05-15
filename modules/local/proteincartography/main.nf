/*
 * modules/local/proteincartography/main.nf
 *
 * Generates an interactive structural landscape map of the viral proteome
 * using ProteinCartography (Arcadia Science).
 * (Bigge et al. Arcadia Science 2024)
 *
 * Runs in "cluster" mode: takes pre-predicted PDB structures, performs
 * all-vs-all Foldseek TM-score comparison, Leiden clustering, and UMAP
 * projection, producing an interactive Plotly HTML map.
 *
 * Pre-filters input PDB directory to pLDDT > params.plddt_cutoff before
 * running (O(N^2) cost makes quality filtering essential at scale).
 *
 * Caveats:
 *   - Best for proteins <1200aa; long/multi-domain structures cluster poorly
 *   - For >1000 structures: ensure >64 GB RAM for all-vs-all TM matrix
 *   - Linux/macOS only (Snakemake + conda envs)
 *
 * Cite: Bigge et al. Arcadia Science 2024
 */

process PROTEINCARTOGRAPHY {
    tag "${meta.id}"
    label 'process_high'

    conda 'conda-forge::python=3.10 bioconda::snakemake=7.32.4 bioconda::foldseek>=9'
    container "${ workflow.containerEngine == 'singularity' ?
        'docker://ghcr.io/ebi-metagenomics/proteincartography:latest' :
        'ghcr.io/ebi-metagenomics/proteincartography:latest' }"

    input:
    tuple val(meta), path(pdb_dir)
    tuple val(meta), path(confidence_tsv)
    tuple val(meta), path(taxonomy_tsv)    // optional overlay; pass empty file if absent

    output:
    tuple val(meta), path("${meta.id}_protein_map/"),             emit: map_dir
    tuple val(meta), path("${meta.id}_protein_map/protein_map.html"), emit: html_map

    when:
    task.ext.when == null || task.ext.when

    script:
    def cutoff  = params.plddt_cutoff ?: 0.7
    def threads = task.cpus
    """
    # Pre-filter PDB dir to pLDDT > cutoff to control O(N^2) cost
    mkdir -p filtered_pdbs

    python3 ${projectDir}/scripts/filter_pdbs_by_plddt.py \\
        --pdb-dir        "${pdb_dir}"      \\
        --confidence-tsv "${confidence_tsv}" \\
        --plddt-cutoff   ${cutoff}         \\
        --outdir         filtered_pdbs

    N=\$(find filtered_pdbs -name "*.pdb" | wc -l)
    echo "Running ProteinCartography on \${N} structures (pLDDT >= ${cutoff})..."

    bash ${projectDir}/scripts/run_proteincartography.sh \\
        filtered_pdbs                    \\
        "${meta.id}_protein_map"         \\
        "${threads}"

    # Overlay taxonomy and pLDDT if taxonomy file is non-empty
    if [ -s "${taxonomy_tsv}" ]; then
        python3 ${projectDir}/scripts/add_taxonomy_overlay.py \\
            --pc-features "${meta.id}_protein_map/features.tsv" \\
            --taxonomy    "${taxonomy_tsv}"                      \\
            --confidence  "${confidence_tsv}"                    \\
            --output      "${meta.id}_protein_map/features_annotated.tsv"
        echo "Taxonomy overlay applied."
    fi
    """

    stub:
    """
    mkdir -p "${meta.id}_protein_map"
    cat > "${meta.id}_protein_map/protein_map.html" << 'HTMLEOF'
<!DOCTYPE html><html><body>
<h2>ProteinCartography stub map — ${meta.id}</h2>
<p>Replace with real ProteinCartography output.</p>
</body></html>
HTMLEOF
    printf "protid\tumap_x\tumap_y\tcluster\n" \
        > "${meta.id}_protein_map/features.tsv"
    """
}
