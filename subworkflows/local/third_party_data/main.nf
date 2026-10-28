nextflow.enable.types = true

// Import modules
include { FALINT                         } from '../../../modules/nf-core/falint/main'
include { PYRODIGAL                      } from '../../../modules/nf-core/pyrodigal/main'
include { GUNZIP as GUNZIP_INPUT         } from '../../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GFF           } from '../../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_FAA           } from '../../../modules/nf-core/gunzip/main'


record ThirdPartyInput {
    meta: Map      // [id, type, source, biome, study_accession, sample_accession]
    fasta: Path
}

record ThirdPartySample {
    meta: Map
    fna: Path
    gff: Path      // predicted by pyrodigal
    faa: Path      // predicted by pyrodigal
}

workflow THIRD_PARTY_DATA {
    take:
    input: Channel<ThirdPartyInput>

    main:
    // Decompress input FASTA files
    compressed = GUNZIP_INPUT(
        input.filter { r -> r.fasta.name.endsWith('.gz') }
            .map { r -> tuple(r.meta, r.fasta) }
    )

    def ch_fasta_ready = input.filter { r -> !r.fasta.name.endsWith('.gz') }
        .map { r -> tuple(r.meta, r.fasta) }
        .mix(compressed.gunzip)

    // FASTA validation
    falint = FALINT(ch_fasta_ready)

    // Keep validated FASTAs
    def ch_output_from_falint = falint.success_log
        .map { meta, _log -> tuple(meta) }
        .join(ch_fasta_ready, by: 0)

    // Report invalid FASTA entries
    invalid_fastas = falint.error_log
        .map { meta, log -> "${meta.id}\t${log}" }
        .collectFile(
            name: "invalid_fastas.tsv",
            newLine: true,
            seed: "sample\tfasta\n",
        )

    // Protein prediction
    pyrodigal = PYRODIGAL(
        ch_output_from_falint,
        'gff'
    )

    faa = GUNZIP_FAA(pyrodigal.faa)
    gff = GUNZIP_GFF(pyrodigal.annotations)

    third_party_data = ch_output_from_falint.map { meta, fasta -> record(meta: meta, fna: fasta) }
        .join(gff.gunzip.map { meta, f -> record(meta: meta, gff: f) }, by: 'meta')
        .join(faa.gunzip.map { meta, f -> record(meta: meta, faa: f) }, by: 'meta')
        .map { r -> record(meta: r.meta, fna: r.fna, gff: r.gff, faa: r.faa) }

    emit:
    third_party_data: Channel<ThirdPartySample> = third_party_data
    invalid_fastas: Channel<Path> = invalid_fastas   // one table, only when some FASTA failed validation
}
