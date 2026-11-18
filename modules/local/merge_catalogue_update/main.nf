nextflow.enable.types = true

/*
 * --mode update: merge the files of the existing catalogue with the update's
 * (bin/merge_catalogue_update.py).
 *
 * MERGE_UPDATE_DROP_IDS  IDs whose rows are dropped from the existing files: the
 *                        representatives absorbed by cluster merges and their proteins
 * MERGE_UPDATE           one file: the existing rows minus the dropped ones, then the update's rows
 * MERGE_UPDATE_REPS      representative IDs of the whole catalogue's cluster TSV
*/

process MERGE_UPDATE_DROP_IDS {

    label 'process_single'
    tag "${meta.id}"
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.75':
        'quay.io/biocontainers/biopython:1.75' }"

    input:
    meta: Map
    dropped_reps: Path     // one old representative per line (VJDB_CLUSTERS_TSV)
    existing_gff: Path     // reps GFF of the existing catalogue

    output:
    file("${meta.id}_drop_ids.txt")

    script:
    """
    merge_catalogue_update.py drop-ids \\
        --dropped-reps ${dropped_reps} \\
        --gff ${existing_gff} \\
        -o ${meta.id}_drop_ids.txt
    """

    stub:
    """
    touch ${meta.id}_drop_ids.txt
    """
}

process MERGE_UPDATE {

    label 'process_single'
    tag "${meta.id}:${key}"
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.75':
        'quay.io/biocontainers/biopython:1.75' }"

    input:
    record(
        meta: Map,
        key: String,           // which file (e.g. 'genomad'), to put the result back in place
        format: String,        // tsv, csv, gff, fasta or tbl
        id_column: String?,    // tsv/csv: ID column name ...
        id_index: Integer?,    // ... or 0-based index
        existing: Path,
        update: Path,
        name: String           // name of the merged file
    )
    drop_ids: Path

    stage:
    stageAs existing, 'existing/*'
    stageAs update, 'update/*'

    output:
    record(
        meta: meta,
        key: key,
        merged: file(name)
    )

    script:
    def id_arg = id_column != null ? "--id-column '${id_column}'" : (id_index != null ? "--id-index ${id_index}" : '')
    """
    merge_catalogue_update.py merge \\
        --format ${format} \\
        ${id_arg} \\
        --old ${existing} \\
        --new ${update} \\
        --drop-ids ${drop_ids} \\
        -o ${name}
    """

    stub:
    """
    touch ${name}
    """
}

process MERGE_UPDATE_REPS {

    label 'process_single'
    tag "${meta.id}"
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.75':
        'quay.io/biocontainers/biopython:1.75' }"

    input:
    meta: Map
    clusters_tsv: Path

    output:
    file("${meta.id}_reps.tsv")

    script:
    """
    merge_catalogue_update.py reps \\
        --clusters-tsv ${clusters_tsv} \\
        -o ${meta.id}_reps.tsv
    """

    stub:
    """
    touch ${meta.id}_reps.tsv
    """
}

/*
 * --mode update: the viral sequence files of the whole catalogue, the existing catalogue's
 * then the update's (viral sequences and prophages: FNA, GFF, FAA). Sequence IDs must be
 * unique; protein IDs are the source IDs and are not checked.
*/
process CONCAT_CATALOGUE_SEQUENCES {

    label 'process_single'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.75':
        'quay.io/biocontainers/biopython:1.75' }"

    input:
    update: Record      // virus_fna, virus_gff, virus_faa, prophage_fna, prophage_gff, prophage_faa
    existing: Record    // the same files of the existing catalogue

    stage:
    stageAs update.virus_fna, 'update/*'
    stageAs update.virus_gff, 'update/*'
    stageAs update.virus_faa, 'update/*'
    stageAs update.prophage_fna, 'update/*'
    stageAs update.prophage_gff, 'update/*'
    stageAs update.prophage_faa, 'update/*'
    stageAs existing.virus_fna, 'existing/*'
    stageAs existing.virus_gff, 'existing/*'
    stageAs existing.virus_faa, 'existing/*'
    stageAs existing.prophage_fna, 'existing/*'
    stageAs existing.prophage_gff, 'existing/*'
    stageAs existing.prophage_faa, 'existing/*'

    output:
    record(
        virus_fna: file('catalogue/viral_sequences.fna'),
        virus_gff: file('catalogue/viral_sequences.gff'),
        virus_faa: file('catalogue/viral_sequences.faa'),
        prophage_fna: file('catalogue/prophages.fna'),
        prophage_gff: file('catalogue/prophages.gff'),
        prophage_faa: file('catalogue/prophages.faa')
    )

    script:
    """
    mkdir catalogue
    merge_catalogue_update.py concat --format fasta --unique-ids -o catalogue/viral_sequences.fna ${existing.virus_fna} ${update.virus_fna}
    merge_catalogue_update.py concat --format gff -o catalogue/viral_sequences.gff ${existing.virus_gff} ${update.virus_gff}
    merge_catalogue_update.py concat --format fasta -o catalogue/viral_sequences.faa ${existing.virus_faa} ${update.virus_faa}
    merge_catalogue_update.py concat --format fasta --unique-ids -o catalogue/prophages.fna ${existing.prophage_fna} ${update.prophage_fna}
    merge_catalogue_update.py concat --format gff -o catalogue/prophages.gff ${existing.prophage_gff} ${update.prophage_gff}
    merge_catalogue_update.py concat --format fasta -o catalogue/prophages.faa ${existing.prophage_faa} ${update.prophage_faa}
    """

    stub:
    """
    mkdir catalogue
    # not empty: empty sequence files are filtered out downstream
    for f in viral_sequences prophages; do
        printf ">stub\\nACGT\\n" > catalogue/\${f}.fna
        printf "##gff-version 3\\n" > catalogue/\${f}.gff
        printf ">stub_1\\nM\\n" > catalogue/\${f}.faa
    done
    """
}
