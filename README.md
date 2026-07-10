# metaViraVerse

<img align="right" width="120" height="120" src="assets/logo.png">

[MGnify](https://www.ebi.ac.uk/metagenomics) Nextflow pipeline to generate **viral catalogue**.

<p align="center">
    <img src="assets/schema.png" alt="Pipeline overview" width="90%">
</p>

## Usage

> [!NOTE]
> This pipeline was written to generate viral catalogue from [MGnify](https://www.ebi.ac.uk/metagenomics/) annotations based on results provided by [emg-viral-pipeline](https://github.com/EBI-Metagenomics/emg-viral-pipeline) (VIRify) and [mobilome-annotation-pipeline](https://github.com/EBI-Metagenomics/mobilome-annotation-pipeline) (MAP).
> 
> If you do not have MGnify results you can still run pipeline using `--third_party_input`.

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
id,gff,fna,faa,type,biome
unique_identifier,viral.gff,viral.fna,viral.faa,metagenome/genome,biome
```

`id` (mandatory) - unique identifier (It is recommended to use ERZ accession if your MAG or assembly was taken ENA) \
`gff` (mandatory) - GFF file containing records in types: _viral_sequence_, _plasmid_, _prophage_. It might also contain CDS records for chosen regions \
`fna` (mandatory) - FASTA file with nucleotide sequences corresponding to chosen regions from GFF \
`faa` (optional) - FASTA file with protein sequences corresponding to CDS regions from GFF \
`type` (mandatory) - string value _genome_ or _metagenome_ describing initial sequence \
`biome` (optional) - metadata describing environmental area of sequence (for example, marine, soil)

## Run

```bash
nextflow run EBI-Metagenomics/metaviraverse \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```


## Citations

If you use this pipeline please make sure to cite all used software.
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
