// Import modules
include { FALINT                         } from '../../../modules/nf-core/falint/main'
include { PYRODIGAL as PYRODIGAL_VIRUS   } from '../../../modules/nf-core/pyrodigal/main'
include { PYRODIGAL as PYRODIGAL_PLASMID } from '../../../modules/nf-core/pyrodigal/main'
include { GUNZIP as GUNZIP_INPUT         } from '../../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GFF_VIRUS     } from '../../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GFF_PLASMID   } from '../../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_FAA_VIRUS     } from '../../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_FAA_PLASMID   } from '../../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_FNA_VIRUS     } from '../../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_FNA_PLASMID   } from '../../../modules/nf-core/gunzip/main'


workflow THIRD_PARTY_DATA {
    take:
    ch_samplesheet // channel: [ val(meta), path(fasta), val(type), val(biome) ]

    main:
    // Decompress input FASTA files
    ch_fasta = ch_samplesheet.map { meta, fasta, _type, _biome -> [meta, fasta] }

    ch_fasta_branched = ch_fasta.branch { _meta, fasta ->
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
        .join(ch_samplesheet.map { meta, _fasta, type, biome -> [meta, type, biome] }, by: 0)
        .map { meta, fasta, type, biome -> [meta, fasta, type, biome ?: 'unknown'] }

    // Report invalid FASTA entries
    FALINT.out.error_log
        .map { meta, log -> "${meta.id}\t${log}" }
        .collectFile(
            name: "invalid_fastas.tsv",
            storeDir: "${params.outdir}",
            newLine: true,
            seed: "sample\tfasta\n",
        )

    // Separate viruses and plasmids
    ch_falint_for_pyrodigal_virus = ch_output_from_falint.filter { _meta, _fasta, type, _biome -> type == 'virus' }
    ch_falint_for_pyrodigal_plasmid = ch_output_from_falint.filter { _meta, _fasta, type, _biome -> type == 'plasmid' }

    // Protein prediction
    PYRODIGAL_VIRUS(
        ch_falint_for_pyrodigal_virus.map { meta, fasta, _type, _biome -> [meta, fasta] },
        'gff',
    )
    PYRODIGAL_PLASMID(
        ch_falint_for_pyrodigal_plasmid.map { meta, fasta, _type, _biome -> [meta, fasta] },
        'gff',
    )

    // Decompress Pyrodigal files
    GUNZIP_FNA_VIRUS(PYRODIGAL_VIRUS.out.fna)
    GUNZIP_FNA_PLASMID(PYRODIGAL_PLASMID.out.fna)
    GUNZIP_GFF_VIRUS(PYRODIGAL_VIRUS.out.annotations)
    GUNZIP_GFF_PLASMID(PYRODIGAL_PLASMID.out.annotations)
    GUNZIP_FAA_VIRUS(PYRODIGAL_VIRUS.out.faa)
    GUNZIP_FAA_PLASMID(PYRODIGAL_PLASMID.out.faa)


    // Re-attach metadata after unzipping
    ch_virus_fna = GUNZIP_FNA_VIRUS.out.gunzip
        .join(ch_falint_for_pyrodigal_virus, by: 0)
        .map { meta, fna, _fasta, type, biome -> [meta, fna, type, biome] }
    ch_plasmid_fna = GUNZIP_FNA_PLASMID.out.gunzip
        .join(ch_falint_for_pyrodigal_plasmid, by: 0)
        .map { meta, fna, _fasta, type, biome -> [meta, fna, type, biome] }

    ch_dedup = ch_virus_fna
        .map { meta, fna, _fasta, _type, biome -> [meta, fna, 'third_party_virus', biome] }
        .mix(ch_plasmid_fna.map { meta, fna, _fasta, _type, biome -> [meta, fna, 'third_party_plasmid', biome] })

    emit:
    fna    = ch_dedup.map { _meta, fna, _type, _biome -> fna }
    faa    = PYRODIGAL_VIRUS.out.faa.mix(PYRODIGAL_PLASMID.out.faa)
    gff    = GUNZIP_GFF_VIRUS.out.gunzip.mix(GUNZIP_GFF_PLASMID.out.gunzip)
    types  = ch_dedup.map { _meta, _fna, type, _biome -> type }
    biomes = ch_dedup.map { _meta, _fna, _type, biome -> biome }
}
