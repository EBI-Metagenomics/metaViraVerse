# metaViraVerse

<img align="right" width="120" height="120" src="assets/logo.png">

A [MGnify](https://www.ebi.ac.uk/metagenomics) Nextflow pipeline for generating a **viral catalogue**.

## Overview

The pipeline takes previously predicted **viral sequences**, **prophages**, and **plasmids** as input and performs downstream analysis and annotation. Viral sequences and prophages are combined into a single **viruses** group, while plasmids are processed as a separate **plasmids** group.

**Viruses:**
- Quality assessment
- rRNA detection
- Clustering at 95% ANI identity
- Taxonomy assignment
- Host detection
- Lifestyle categorisation
- Protein prediction
- Protein annotation
- Protein clustering

**Plasmids:**
- Clustering at 85% ANI identity
- Host detection [in development]
- Prediction of replicon family, relaxase type, and mate-pair formation type [in development]

<p align="center">
    <img src="assets/schema.png" alt="Pipeline overview" width="90%">
</p>

## Input

> [!NOTE]
> This pipeline was originally written to build viral catalogues from [MGnify](https://www.ebi.ac.uk/metagenomics/) annotations, using results produced by [emg-viral-pipeline](https://github.com/EBI-Metagenomics/emg-viral-pipeline) (VIRify) and [mobilome-annotation-pipeline](https://github.com/EBI-Metagenomics/mobilome-annotation-pipeline) (MAP).
>
> If you don't have MGnify results, you can still run the pipeline using `--third_party_input`.

First, prepare a samplesheet describing your input data:

`samplesheet.csv`:

```csv
id,gff,fna,faa,type,biome
unique_identifier,viral.gff,viral.fna,viral.faa,metagenome/genome,biome
```

| Column  | Required | Description |
|---------|----------|-------------|
| `id`    | Yes | Unique identifier. We recommend using the ERZ accession if the MAG or assembly originates from ENA. |
| `gff`   | Yes | GFF file containing viral and plasmid records. May also contain CDS records for the selected regions. |
| `fna`   | Yes | FASTA file with nucleotide sequences for the selected regions in the GFF. |
| `faa`   | No | FASTA file with protein sequences for the CDS regions in the GFF. |
| `type`  | Yes | `genome` (for a MAG source) or `metagenome` (for an assembly source), describing the origin of the sequence. |
| `biome` | No | Metadata describing the sequence's environmental origin (for example: marine, soil). |

## Usage

- For MGnify input, see the [MGnify usage guide](docs/mgnify_usage.md).
- For third-party data, see the [third-party usage guide](docs/third_party_usage.md).

## Run

```bash
nextflow run EBI-Metagenomics/metaviraverse \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

## Citations

If you use this pipeline, please cite all software it uses.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **MGnify: the microbiome sequence data analysis resource in 2023**
>
> Richardson L, Allen B, Baldi G, Beracochea M, Bileschi ML, Burdett T, et al.
>
> Vol. 51, Nucleic Acids Research. Oxford University Press (OUP); 2022. p. D753–9. Available from: http://dx.doi.org/10.1093/nar/gkac1080
>
> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
