# metaViraVerse

<img align="right" width="120" height="120" src="assets/logo.png">

[MGnify](https://www.ebi.ac.uk/metagenomics) Nextflow pipeline to generate **viral catalogue** from assemblies.

<p align="center">
    <img src="assets/schema.png" alt="Pipeline overview" width="90%">
</p>

## Usage

> [!NOTE]
> This pipeline is based on results provided by [emg-viral-pipeline](https://github.com/EBI-Metagenomics/emg-viral-pipeline) (VIRify). In order to generate a catalogue you need to launch VIRify on each sequence file in advance.

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
id,gff,fasta
unique_identifier,assembly_virify.gff,assembly_virify.fasta
```

Each row represents an sequence file that was annotated with [VIRify](https://github.com/EBI-Metagenomics/emg-viral-pipeline).

## Run

```bash
nextflow run EBI-Metagenomics/metaviraverse \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

## Output

```
├── plasmids
├── viral_sequences
├── pipeline_info
```

Folder `plasmids`:

```
├── cluster_reps
 ──── plasmids_reps.fasta.gz    # cluster rep compressed fasta

├── clustering
 ──── plasmids_clusters.tsv     # clusters (representative \t members)
 ──── plasmids_pairani.tsv      # blastn pairani table

├── plasmids.fasta              # all sequences
```

Folder `viral_sequences`:

```
├── cluster_reps

 ──── crisprcasfinder
 ──────── viral_sequences_crisprcasfinder.gff
 ──────── viral_sequences_crisprcasfinder.tsv
 ──────── viral_sequences_crisprcasfinder_hq.gff

 ──── taxonomy_plot                          # ViPhOG taxonomy
 ──────── viral_sequences_krona.html         # krona plot
 ──────── viral_sequences_krona.tsv          # taxonomy table with counts (count \t taxonomy tav-separated)
 ──────── viral_sequences_sankey.html        # sankey plot

 ──── taxonomy_vitap                         # ICTV taxonomy
 ──────── viral_sequences_vitap_best.tsv     # taxonomy table with counts (count \t taxonomy tav-separated)
 ──────── viral_sequences_vitap_sankey.html  # sankey plot

 ──── viral_sequences_reps.fasta.gz          # cluster rep compressed fasta
 ──── viral_sequences_reps_stats.tsv         # basic statistics calculated per each cluster representative

├── clustering
 ──── viral_sequences_clusters.tsv           # clusters (representative \t members)
 ──── viral_sequences_pairani.tsv            # blastn pairani table

├── viral_sequences.fasta                    # all sequences
```

### Under review:

**VITAP** (https://www.nature.com/articles/s41467-025-57500-7)
Database: https://figshare.com/articles/dataset/The_database_of_VITAP_2024_03_18_/25426159/3?file=49682337
taxonomy + sankey plot

**PLSDB** (https://academic.oup.com/nar/article/53/D1/D189/7905312)
plasmids only screening with mash

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
