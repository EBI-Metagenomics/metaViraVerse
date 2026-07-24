// Import modules
include { FALINT                         } from '../../../modules/nf-core/falint/main'
include { PYRODIGAL                      } from '../../../modules/nf-core/pyrodigal/main'
include { GUNZIP as GUNZIP_INPUT         } from '../../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GFF_VIRUS     } from '../../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GFF_PLASMID   } from '../../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_FAA_VIRUS     } from '../../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_FAA_PLASMID   } from '../../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_FNA_VIRUS     } from '../../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_FNA_PLASMID   } from '../../../modules/nf-core/gunzip/main'
include { PLASQUID_WORKFLOW                       } from '../plasquid_workflow'


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
    plasmid_proteins = PYRODIGAL.out.faa.filter{meta, faa -> meta.type == 'plasmid'}
    plasmid_fna      = PYRODIGAL.out.fna.filter{meta, fna -> meta.type == 'plasmid'}

    GUNZIP_FAA_PLASMID(plasmid_proteins)
    GUNZIP_FNA_PLASMID(plasmid_fna)

    // TODO remove after testing
    PLASQUID_WORKFLOW(
        GUNZIP_FNA_PLASMID.out.gunzip,
        GUNZIP_FAA_PLASMID.out.gunzip
    )
}