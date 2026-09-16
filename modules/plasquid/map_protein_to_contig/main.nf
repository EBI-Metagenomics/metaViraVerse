process MAP_PROTEIN_TO_CONTIG {
    /*
     * Every plasquid script needs to resolve a protein/RNA hit back to its parent
     * contig. For RNA hits (cmsearch, against the nucleotide fasta) the hit id
     * already *is* the contig id. For protein hits it isn't: this pipeline's
     * protein FASTA keeps each ORF's original (pre-rename) id, e.g.
     * "MGYG000517684_26|plasmid-1:6000_1", while the matching nucleotide contig has
     * since been renamed to a short accession, e.g. "seq15" -- the two no longer
     * share a common prefix, so it can't be recovered by string-stripping the
     * protein id alone.
     *
     * Build the crosswalk directly from the representative GFF instead: every CDS
     * feature's own seqid (column 1, already renamed) paired with its `ID=`
     * attribute (the untouched protein id) gives an exact (protein_id, contig) map.
    */

    label 'process_single'
    tag "${meta.id}"
    container "docker://mgimenez720/plasquid:latest"

    input:
    tuple val(meta), path(gff)

    output:
    tuple val(meta), path("protein_to_contig.tsv"), emit: map
    path "versions.yml",                             emit: versions

    script:
    """
    printf "protein_id\\tcontig\\n" > protein_to_contig.tsv
    awk -F'\\t' '\$3 == "CDS" { match(\$9, /ID=[^;]+/); id = substr(\$9, RSTART + 3, RLENGTH - 3); print id"\\t"\$1 }' ${gff} >> protein_to_contig.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        protein_to_contig_map: "gff-derived (awk)"
    END_VERSIONS
    """
}
