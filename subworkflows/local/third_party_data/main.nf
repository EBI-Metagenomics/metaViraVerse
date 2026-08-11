// Import modules
include { FALINT                         } from '../../../modules/nf-core/falint/main'
include { PYRODIGAL                      } from '../../../modules/nf-core/pyrodigal/main'
include { GUNZIP as GUNZIP_INPUT         } from '../../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GFF           } from '../../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_FAA           } from '../../../modules/nf-core/gunzip/main'


workflow THIRD_PARTY_DATA {
    take:
    ch_input // channel: [ val(meta), path(fasta) ]

    main:
    // Decompress input FASTA files

    ch_fasta_branched = ch_input.branch { _meta, fasta ->
        compressed: fasta.name.endsWith('.gz')
        uncompressed: true
    }

    GUNZIP_INPUT(ch_fasta_branched.compressed)

    ch_fasta_ready = ch_fasta_branched.uncompressed.mix(GUNZIP_INPUT.out.gunzip)

    // FASTA validation
    FALINT(ch_fasta_ready)

    // Keep validated FASTAs
    ch_output_from_falint = FALINT.out.success_log
        .map { meta, _log -> meta }
        .join(ch_fasta_ready, by: 0)

    // Report invalid FASTA entries
    FALINT.out.error_log
        .map { meta, log -> "${meta.id}\t${log}" }
        .collectFile(
            name: "invalid_fastas.tsv",
            storeDir: "${params.outdir}",
            newLine: true,
            seed: "sample\tfasta\n",
        )

    // Protein prediction
    PYRODIGAL(
        ch_output_from_falint,
        'gff'
    )

    GUNZIP_FAA(PYRODIGAL.out.faa)
    GUNZIP_GFF(PYRODIGAL.out.annotations)

    emit:
    third_party_data = ch_output_from_falint.join(GUNZIP_GFF.out.gunzip).join(GUNZIP_FAA.out.gunzip)  // [meta, fasta, gff, faa]
}
