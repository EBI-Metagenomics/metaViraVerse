// Import modules
include { FALINT                         } from '../../../modules/nf-core/falint/main'
include { PYRODIGAL as PYRODIGAL_VIRUS   } from '../../../modules/nf-core/pyrodigal/main'
include { PYRODIGAL as PYRODIGAL_PLASMID } from '../../../modules/nf-core/pyrodigal/main'
include { GUNZIP as GUNZIP_VIRUS         } from '../../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_PLASMID       } from '../../../modules/nf-core/gunzip/main'

workflow THIRD_PARTY_DATA {
    take:
    ch_samplesheet // channel: [ val(meta), path(fasta), val(type), val(biome) ]

    main:

    // FASTA validation
    FALINT(
        ch_samplesheet.map { meta, fasta, type, biome -> [meta, fasta] }
    )

    // Keep validated FASTAs
    ch_validated = FALINT.out.success_log
        .join(ch_samplesheet, by: 0)
        .map { success_log, meta, fasta, type, biome -> [meta, fasta, type, biome ?: 'unknown'] }

    // Separate viruses and plasmids
    ch_viruses = ch_validated.filter { meta, fasta, type, biome -> type == 'virus' }
    ch_plasmids = ch_validated.filter { meta, fasta, type, biome -> type == 'plasmid' }

    // Protein prediction
    PYRODIGAL_VIRUS(
        ch_viruses.map { meta, fasta, type, biome -> [meta, fasta] },
        'fna',
    )
    PYRODIGAL_PLASMID(
        ch_plasmids.map { meta, fasta, type, biome -> [meta, fasta] },
        'fna',
    )
    
    // Decompress FASTA files
    GUNZIP_VIRUS(
            PYRODIGAL_VIRUS.out.fna
        )
    GUNZIP_PLASMID(
        PYRODIGAL_PLASMID.out.fna
    )


    // Re-attach metadata after unzipping
    ch_virus_fna = GUNZIP_VIRUS.out.gunzip
        .join(ch_viruses, by: 0)
        .map { meta, fna, fasta, type, biome -> [meta, fna, type, biome] }
    ch_plasmid_fna = GUNZIP_PLASMID.out.gunzip
        .join(ch_plasmids, by: 0)
        .map { meta, fna, fasta, type, biome -> [meta, fna, type, biome] }

    ch_dedup = ch_virus_fna.mix(ch_plasmid_fna)

    emit:
    fna      = ch_dedup.map { meta, fna, type, biome -> fna }
    faa      = PYRODIGAL_VIRUS.out.faa.mix(PYRODIGAL_PLASMID.out.faa)
    types    = ch_dedup.map { meta, fna, type, biome -> type }
    biomes   = ch_dedup.map { meta, fna, type, biome -> biome }
}
