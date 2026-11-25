# Questions and answers

Common questions about the viral and plasmid catalogues produced by metaViraVerse. For the description of every file and metadata column see [Output](output.md); for how each step works see [Methods](methods.md); for tool and database versions see [Databases](databases.md).

If your question is not answered here, please open an issue at <https://github.com/EBI-Metagenomics/metaViraVerse/issues>.

## Contents

- [Identifiers](#identifiers)
- [Sequences and clusters](#sequences-and-clusters)
- [Missing values](#missing-values)
- [Taxonomy and hosts](#taxonomy-and-hosts)
- [Plasmids](#plasmids)

## Identifiers

### Why were all input sequences renamed?

Every sequence in the catalogue gets a short, stable identifier such as `MGYV0000000123`. Sequences are collected from many MAG catalogues, assembly analyses and, optionally, third-party sources, and their original names are long, follow different conventions and are not guaranteed to be unique across sources. Several of the tools used by the pipeline also fail on long or unusual FASTA headers. A single numbering scheme gives every sequence one unique identifier that can be used in the MGnify database and across all files of the catalogue. See [Rename sequences](methods.md#rename-sequences).

### How do I find the original name of a sequence?

Look up its identifier (`ID`) in `viruses/viruses-all-metadata.tsv.gz` or `plasmids/plasmids-all-metadata.tsv.gz`: the `Original_ID` column has the original name. `additional_data/rename_mapfile.tsv` maps every identifier to its original name too, including sequences that were later excluded.

### How do I find the genome or assembly a sequence was found in?

The original name starts with the accession of its source, followed by the contig number, the sequence type and its coordinates on the contig. For example, `MGYG000519004_12|viral_sequence-1:32972` is a viral sequence at positions 1–32,972 of contig 12 of genome `MGYG000519004`.

- `MGYG...` is an MGnify genome (MAG or isolate). Search for it on the [MGnify genomes](https://www.ebi.ac.uk/metagenomics/browse/genomes) website. The `Source_accession` column of the metadata tables gives its INSDC genome accession (e.g. `GCA_...`), and `Source_sample` and `Source_project` give its sample and study.
- `ERZ...` is an assembly analysed by MGnify. Search for it on [MGnify](https://www.ebi.ac.uk/metagenomics) or the [ENA browser](https://www.ebi.ac.uk/ena/browser/). Its sample and study are in the `Source_sample` and `Source_project` columns.

## Sequences and clusters

### Why is the number of sequences in the catalogue lower than the number of input sequences?

Two steps remove sequences before clustering (see [Sequences deduplication](methods.md#sequences-deduplication)):

- **Duplicates.** Identical sequences, for example the same virus found in an assembly and in a MAG built from it, are kept only once. The biomes of all copies are merged into the kept record.
- **Quality control.** Viral sequences and prophages are removed when they carry an rRNA gene (a sign of host contamination), or when CheckV could not determine their quality, found no viral genes and reports a high k-mer frequency. Plasmids are not filtered.

Every removed sequence is listed in `additional_data/catalogue_excluded_records_metadata.tsv`, with the reason in the `filter_reason` column (`duplicate`, `rrna` or `not_determined`) and its original name in `original_name`. Some of the sequences removed by quality control may still be novel viruses, so the table can be used to go back to the source genome or assembly and inspect them. The numbers are in `catalogue.json`: `total_sequences` (input), `duplicate_sequences`, `rrna_excluded_sequences`, `not_determined_excluded_sequences` and `unique_sequences` (kept). Deduplication runs on all categories pooled together, and the kept copy decides the category, so the drop within one category (e.g. plasmids, which are never removed by quality control) can come from duplicates alone.

### Why are there sequences with CheckV quality `Not-determined` in the catalogue?

CheckV quality `Not-determined` alone is not a reason to remove a sequence. CheckV reports it when it can't estimate completeness, which is common for short fragments and for novel viruses with no close reference. A viral sequence or prophage is removed as `not_determined` only when all three of these hold:

1. `checkv_quality` is `Not-determined`;
2. CheckV found no viral genes (`viral_genes` = 0);
3. the k-mer frequency is 1.0 or higher (`kmer_freq` ≥ 1.0).

No viral genes together with a high k-mer frequency is typical of host contamination. A `Not-determined` sequence that has at least one viral gene is kept as a real virus that CheckV could not assess. Plasmids are never filtered by CheckV quality.

The final metadata table has `checkv_quality` and `kmer_freq`, but not `viral_genes`. To see why a sequence was kept, look up its `viral_genes` in the metadata table in `additional_data/choose_sequences/`, which keeps all CheckV columns. Kept `Not-determined` sequences usually have low or `NA` completeness, so filter on `completeness` or `miuvig_quality` if you need only well-characterised genomes.

### Which sequences are annotated?

Only cluster representatives: one sequence per viral cluster and one per plasmid cluster. Viral sequences and prophages are clustered at 95% average nucleotide identity (ANI) and 85% coverage, roughly species level. Plasmids are clustered separately. Taxonomy, host, lifestyle and protein annotations are therefore filled in for representatives only. To get them for any other sequence, look up its `Cluster_rep` in the metadata table and use the values of that representative. When a catalogue is updated with new sequences, only the representatives of new clusters are annotated.

### Why does the GFF contain more annotations than the tools run by the pipeline?

The input GFF files come from the MGnify genome and assembly analysis pipelines, which had already annotated the proteins (for example with InterPro, Pfam or eggNOG). These annotations are kept on the CDS records. metaViraVerse adds only virus- and plasmid-specific annotations: viral taxonomy, hosts, lifestyle, HMMER searches against viral and other HMM databases, antimicrobial resistance genes, and plasmid typing and mobility. The full list is in the README of the catalogue and in [Methods](methods.md).

## Missing values

### Why do some metadata columns contain `missing`, `NA` or `not-provided`?

These values mean the information is not available, for one of these reasons:

- **The step was skipped.** For example, rRNA is `not-provided` when rRNA detection was skipped, and the CheckV columns are `NA` when CheckV was skipped.
- **The sequence is not a cluster representative.** Taxonomy and host columns are filled for representatives only (see [Which sequences are annotated?](#which-sequences-are-annotated)).
- **The tool gave no result.** For example, a sequence without a known relative may have no taxonomic assignment, and iPHoP may find no host.
- **The source has no such metadata.** For example, sequences from assemblies have no source genome lineage (see [Why is the source lineage missing for some sequences?](#why-is-the-source-lineage-missing-for-some-sequences)), and third-party sequences may come without a sample or study.

`missing` is written by the final metadata step. `NA` and `not-provided` come from earlier steps that already had no value. For every column, [Output](output.md#metadata) says whether it can be missing and when.

## Taxonomy and hosts

### Why do geNomad, VITAP and ViPhOGs give different lineages for the same sequence?

The three tools classify sequences in different ways and use different reference data:

- **geNomad** uses marker genes and a neural-network model trained on reference viral genomes.
- **VITAP** compares the sequence with reference viral genomes at the nucleotide and protein level.
- **ViPhOGs** uses the pattern of hits to viral protein HMMs along the genome.

Each tool reports a lineage only down to the rank it is confident about, so one tool may stop at class while another reaches genus. Their reference databases also follow different releases of the viral taxonomy. The pipeline reports all three side by side and does not combine them into a consensus, so you can compare them or choose the one that best fits your question. Agreement between tools increases confidence in a lineage. See [Taxonomy assignment (DNA sequence based)](methods.md#taxonomy-assignment-dna-sequence-based) and [Taxonomy assignment (Aminoacid sequence based)](methods.md#taxonomy-assignment-aminoacid-sequence-based).

### How reliable are the host predictions?

They are predictions, and they are only as good as the reference data behind them:

- **iPHoP** combines several kinds of host evidence (CRISPR spacer matches, sequence similarity to host genomes, k-mer composition) and reports only predictions that pass its confidence cut-off. Where it lists several possible hosts, the metadata keeps the one with the highest confidence score. Hosts that are not in the iPHoP database cannot be predicted; the database version is in [Databases](databases.md).
- **SpacePHARER** (when used) matches viral sequences to CRISPR spacers from the genomes in the MGnify catalogues. It can only find hosts among those genomes. A match is strong evidence of a past infection of that host lineage. Where there are several matches, the one with the lowest e-value is kept.

A missing host means no prediction passed these criteria, not that the virus has no host. Genus-level predictions are generally more reliable than genome-level ones.

### Why is the source lineage missing for some sequences?

`Source_lineage` is the taxonomy of the genome (MAG) the sequence was found in. Sequences found in assemblies come from a mix of organisms, not from a single genome, so they have no source lineage. The same applies to third-party sequences without source information. For viruses, the source lineage is the lineage of the genome carrying the sequence, which for a prophage is usually its host. It is not a host prediction for free viral sequences; use the host columns for that.

## Plasmids

### What is the host of a plasmid?

The pipeline does not predict plasmid hosts. The genome (MAG) a plasmid was found in is treated as its host, so the host lineage is the `Source_lineage` column of `plasmids/plasmids-all-metadata.tsv.gz`. Plasmids found in assemblies have no source genome, so their host is `missing` (see [Why is the source lineage missing for some sequences?](#why-is-the-source-lineage-missing-for-some-sequences)). MOB-suite also predicts a likely host range from the plasmid's replicon and relaxase types. It is in `mobsuite-typer/plasmids_report.txt` and not in the metadata table.
