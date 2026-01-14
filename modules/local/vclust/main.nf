process VCLUST {

    label 'process_low'
    tag "${meta.id}"
    container 'community.wave.seqera.io/library/pip_vclust:0dc3bb66940e8323'

    input:
      tuple val(meta), path(fasta)
      val(ani_threshold)            // 0 < value < 1
      val(coverage_threshold)       // 0 < value < 1

    output:
      tuple val(meta), path("${meta.id}_clusters.tsv"),      emit: clusters_tsv
      tuple val(meta), path("${meta.id}_fltr.txt"),          emit: filter_txt
      tuple val(meta), path("${meta.id}_ani.tsv"),           emit: align_ani_tsv
      tuple val(meta), path("${meta.id}_ani.ids.tsv"),       emit: align_ani_ids
      path "versions.yml",                                   emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    # Deduplicate
    echo "Run deduplicate"
    vclust deduplicate -i ${fasta} -o ${meta.id}_dedup.fasta

    # Prefilter similar genome sequence pairs before conducting pairwise alignments.
    # Create pre-alignment filter with 20 common 25-mers and 70% identity over the shorter sequence.
    # Process genomes in batches of 2 million sequences
    # TODO: if needed: Limit the number of target sequences to top 1000 per query genome with --max-seqs
    echo "Run prefilter"
    vclust prefilter \
      -i ${meta.id}_dedup.fasta \
      -o ${meta.id}_fltr.txt \
      --min-kmers 20 \
      --min-ident ${ani_threshold} \
      --batch-size 2000000

    # Align similar genome sequence pairs and calculate pairwise ANI measures.
    # ANI ≥ ani_threshold and query coverage ≥ coverage_threshold
    # To manage the output TSV file size, Vclust offers three formats: standard, lite, complete
    echo "Run align"
    vclust align \
      -i ${meta.id}_dedup.fasta \
      -o ${meta.id}_ani.tsv \
      --filter ${meta.id}_fltr.txt \
      --out-ani ${ani_threshold} \
      --out-qcov ${coverage_threshold} \
      --outfmt standard

    # Cluster genome sequences based on given ANI measure and minimum threshold.
    echo "Run cluster"
    vclust cluster \
      -i ${meta.id}_ani.tsv \
      -o ${meta.id}_clusters.tsv \
      --ids ${meta.id}_ani.ids.tsv \
      --metric ani \
      --ani ${ani_threshold} \
      --qcov ${coverage_threshold} \
      --rcov ${coverage_threshold} \
      --algorithm leiden \
      --out-repr

    echo "Done."

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vclust: \$(vclust --version 2>&1 | sed 's/v //g')
    END_VERSIONS
    """
}
