# EBI-Metagenomics/metaviraverse: Output

## Introduction

This document describes the output produced by the pipeline. The pipeline builds a viral catalogue and a plasmid catalogue from predicted viral sequences, prophages and plasmids. Results are organised under two top-level directories, **viruses/** and **plasmids/**, plus supporting directories for intermediate data, QC and pipeline metadata.

The directories listed below are created in the results directory when the pipeline finishes. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline produces the following outputs:

- [Metadata](#metadata) - combined metadata for the viral and plasmid catalogues
- [Viruses](#viruses) - viral sequences and prophages, clustered and annotated with taxonomy, host, lifestyle and function
- [Plasmids](#plasmids) - plasmid sequences, clustered and annotated with AMR genes and MOB typing
- [Catalogue summary](#catalogue-summary) - top-level statistics of the catalogue
- [Additional data](#additional-data) - intermediate files produced while preparing the input (renaming, quality control, rRNA detection, sequence splitting)
- [MultiQC](#multiqc) - aggregate report describing results and QC from the whole pipeline
- [Pipeline information](#pipeline-information) - reports and metrics generated during the workflow execution

## Metadata

`viruses-all-metadata.tsv.gz` and `plasmids-all-metadata.tsv.gz` list every sequence in the catalogue. Each row contains the sequence identifiers, the sample, project and biome it came from, the taxonomic lineage of the source genome, and the sequence length and checksum. For viruses, the table also includes CheckV quality metrics and the taxonomy assigned by VITAP, ViPhOGs and geNomad. `viruses-cluster-stats.tsv.gz` summarises each viral cluster, one row per cluster representative.

<details markdown="1">
<summary>Output files</summary>

- `viruses/viruses-all-metadata.tsv.gz`: one row per **input** viral sequence or prophage (before clustering), combining source, quality and taxonomy information.
- `viruses/viruses-cluster-stats.tsv.gz`: one row per viral **cluster representative**, with per-cluster statistics (cluster size, mean gene and viral-gene counts, biomes and source types in the cluster, completeness and contamination ranges).
- `plasmids/plasmids-all-metadata.tsv.gz`: one row per **input** plasmid sequence (before clustering).

</details>

Missing values are written as `missing` by the metadata step, or as `NA` / `not-provided` when the value was already absent in an earlier step. In the tables below, **Can be missing** says whether a column can hold one of these values instead of real data, and when.

#### `viruses-all-metadata.tsv.gz`

| Column | Can be missing | Description |
| --- | --- | --- |
| `ID` | No | Catalogue identifier assigned by the pipeline (e.g. `MGYV0000000001`). |
| `Original_ID` | No | Original sequence name from the input, e.g. `MGYG000519004_12\|viral_sequence-1:32972` (source genome, contig, sequence type and coordinates). |
| `Source_of_prediction` | Yes (`NA`) | Tool that predicted the sequence, taken from the source column of its GFF record (e.g. `geNomad`); `third_party` for sequences from `--third_party_input`. |
| `rRNA` | Yes (`not-provided`) | `Yes` if barrnap found an rRNA gene on the sequence, otherwise `No`; `not-provided` when rRNA detection was skipped (`--skip_rrna_detection`). |
| `Cluster_rep` | Yes (`missing`) | `ID` of the representative of the cluster this sequence belongs to; equal to `ID` for representatives. |
| `vitap_lineage` | Yes (`missing`) | Taxonomy assigned by VITAP. Only cluster representatives are classified, so other members are `missing`; also `missing` when VITAP gave no assignment or was skipped. |
| `viphogs_lineage` | Yes (`missing`) | Taxonomy assigned from ViPhOG HMM hits, semicolon-separated from superkingdom to genus. Same rules as `vitap_lineage`. |
| `genomad_lineage` | Yes (`missing`) | Taxonomy assigned by geNomad. Same rules as `vitap_lineage`. |
| `Host_iphop_genome_lineage` | Yes (`missing`) | Host lineage (GTDB, down to species) of the highest-confidence iPHoP genome-level prediction. `missing` when iPHoP found no host or was skipped (`--skip_iphop`). |
| `Host_iphop_genus_lineage` | Yes (`missing`) | Host lineage (GTDB, down to genus) of the highest-confidence iPHoP genus-level prediction. Same rules as `Host_iphop_genome_lineage`. |
| `Host_spacepharer_lineage` | Yes (`missing`) | Host lineage from the best (lowest e-value) SpacePHARER CRISPR spacer match. Only filled when `--predict_host_from_custom_spacers` is used. |
| `Source` | Yes (`NA`) | Type of the source the sequence came from: `genome` (isolate genome) or `metagenome` (MAG). |
| `Biome` | Yes (`NA`) | Biome(s) of the input catalogue(s) the sequence was found in, comma-separated when an identical sequence came from several. |
| `Source_accession` | Yes (`missing`) | Accession of the source genome (e.g. `GCA_001650725`). `missing` for third-party sequences or genomes absent from the catalogue metadata. |
| `Source_lineage` | Yes (`missing`) | GTDB lineage of the source genome (the MAG the sequence was found in). Same rules as `Source_accession`. |
| `Source_sample` | Yes (`missing`) | ENA/BioSample accession of the sample the source genome came from. Same rules as `Source_accession`. |
| `Source_project` | Yes (`missing`) | ENA study/project accession of the source genome. Same rules as `Source_accession`. |
| `Sequence_length` | No | Sequence length in bp. |
| `Sequence_sha256` | No | SHA-256 checksum of the (upper-case) nucleotide sequence, used to detect duplicates. |
| `checkv_quality` | Yes (`NA`) | CheckV quality tier (`Complete`, `High-quality`, `Medium-quality`, `Low-quality`, `Not-determined`). All CheckV columns are `NA` when CheckV was skipped (`--skip_checkv`). |
| `miuvig_quality` | Yes (`NA`) | Quality tier following the MIUViG standard (`High-quality`, `Genome-fragment`, ...). |
| `completeness` | Yes (`NA`) | Estimated completeness (%); `NA` when CheckV could not estimate it. |
| `completeness_method` | Yes (`NA`) | Method CheckV used to estimate completeness (e.g. AAI-based, HMM-based). |
| `contamination` | Yes (`NA`) | Estimated contamination (%): the share of the sequence identified as host-derived. |
| `provirus` | Yes (`NA`) | `Yes` if CheckV identified the sequence as a provirus (with flanking host regions), otherwise `No`. |
| `proviral_length` | Yes (`NA`) | Length (bp) of the viral region when `provirus` is `Yes`. |
| `kmer_freq` | Yes (`NA`) | Average k-mer frequency; values well above 1 indicate a sequence with repeated regions (possible assembly artefact). |

#### `viruses-cluster-stats.tsv.gz`

One row per viral cluster, described by its representative and by aggregates over all cluster members.

| Column | Can be missing | Description |
| --- | --- | --- |
| `Rep_ID` | No | `ID` of the cluster representative. |
| `Rep_vitap_lineage` | Yes (`missing`) | VITAP taxonomy of the representative; `missing` when there was no assignment or VITAP was skipped. |
| `Rep_viphogs_lineage` | Yes (`missing`) | ViPhOGs taxonomy of the representative. Same rules as `Rep_vitap_lineage`. |
| `Rep_genomad_lineage` | Yes (`missing`) | geNomad taxonomy of the representative. Same rules as `Rep_vitap_lineage`. |
| `Rep_source` | Yes (`missing`) | Source type of the representative: `genome` or `metagenome`. |
| `Rep_biome` | Yes (`missing`) | Biome(s) of the representative, comma-separated. |
| `cluster_size` | No | Number of sequences in the cluster, including the representative. |
| `cluster_mean_viral_genes` | Yes (`missing`) | Mean number of viral genes (from CheckV) across cluster members; `missing` when no member has CheckV data. |
| `cluster_mean_gene_count` | Yes (`missing`) | Mean number of genes (from CheckV) across cluster members. Same rules as `cluster_mean_viral_genes`. |
| `cluster_biomes` | Yes (`missing`) | All biomes found across cluster members, semicolon-separated. |
| `cluster_types` | Yes (`missing`) | All source types found across cluster members (`genome`, `metagenome`), semicolon-separated. |
| `cluster_completeness_range` | Yes (`missing`) | Range of CheckV completeness across members (`min-max`, or a single value when only one member has an estimate). |
| `cluster_contamination_range` | Yes (`missing`) | Range of CheckV contamination across members, in the same format as `cluster_completeness_range`. |

#### `plasmids-all-metadata.tsv.gz`

Plasmids are not assessed with CheckV and are not assigned viral taxonomy or hosts, so this table has a subset of the virus columns.

| Column | Can be missing | Description |
| --- | --- | --- |
| `ID` | No | Catalogue identifier assigned by the pipeline. |
| `Original_ID` | No | Original sequence name from the input (source genome, contig, sequence type and coordinates). |
| `Cluster_rep` | Yes (`missing`) | `ID` of the representative of the cluster this plasmid belongs to; equal to `ID` for representatives. |
| `Source` | Yes (`NA`) | Type of the source the plasmid came from: `genome` or `metagenome`. |
| `Biome` | Yes (`NA`) | Biome(s) of the input catalogue(s) the plasmid was found in, comma-separated. |
| `Source_accession` | Yes (`missing`) | Accession of the source genome. `missing` for third-party sequences or genomes absent from the catalogue metadata. |
| `Source_lineage` | Yes (`missing`) | GTDB lineage of the source genome, which is also treated as the plasmid host. Same rules as `Source_accession`. |
| `Source_sample` | Yes (`missing`) | ENA/BioSample accession of the sample the source genome came from. Same rules as `Source_accession`. |
| `Source_project` | Yes (`missing`) | ENA study/project accession of the source genome. Same rules as `Source_accession`. |
| `Sequence_length` | No | Sequence length in bp. |
| `Sequence_sha256` | No | SHA-256 checksum of the (upper-case) nucleotide sequence. |

## Viruses

Viral sequences and prophages predicted upstream (for example by [VIRify](https://github.com/EBI-Metagenomics/emg-viral-pipeline) or the [mobilome-annotation-pipeline](https://github.com/EBI-Metagenomics/mobilome-annotation-pipeline)) are pooled into a single **viruses** set and clustered, and the cluster representatives are annotated.

<details markdown="1">
<summary>Output files</summary>

- `viruses/`
  - `prophages.fna.gz`, `prophages.faa.gz`, `prophages.gff.gz`: nucleotide sequences, proteins and GFF annotation of all prophages.
  - `viral_sequences.fna.gz`, `viral_sequences.faa.gz`, `viral_sequences.gff.gz`: nucleotide sequences, proteins and GFF annotation of all viral sequences.
  - `viruses_host_sankey.html`: Sankey plot of the taxonomy of the source genomes (MAGs) the viruses were found in.
  - `viruses-all-metadata.tsv.gz`: one row per **input** viral sequence or prophage (before clustering), combining source, quality and taxonomy information.
  - `viruses-cluster-stats.tsv.gz`: one row per viral **cluster representative**, with per-cluster statistics (cluster size, mean gene and viral-gene counts, biomes and source types in the cluster, completeness and contamination ranges).
  - `cluster_representatives/`: results for the cluster representatives. Sequences are clustered with [vclust](https://github.com/refresh-bio/vclust) at 95% ANI and 85% coverage.
    - `viruses_reps.fasta.gz`, `viruses_reps.faa.gz`: nucleotide and protein sequences of one representative per viral cluster.
    - `viruses_reps_annotated.gff.gz`: GFF of the representatives, enriched with lifestyle (BACPHLIP), HMMER (e.g. PVOG, TIGRFAM) and AMR annotations.
    - `functional_annotation/`
      - `amr/`
        - `viruses.tsv`: antimicrobial resistance calls from [AMRFinderPlus](https://github.com/ncbi/amr), [DeepARG](https://github.com/gaarangoa/deeparg) and [RGI](https://github.com/arpcard/rgi).
        - `integrated_viruses.gff`: GFF with the AMR annotations merged in.
      - `hmmer/`: [HMMER](http://hmmer.org/) hits against HMM profile databases (`*.tbl.gz`), with a per-accession function summary (`*_summary.tsv.gz`).
      - `lifestyle/`: virulent/temperate lifestyle prediction from [BACPHLIP](https://github.com/adamhockenberry/bacphlip) (`viruses.bacphlip.gz`).
    - `host_detection/`
      - `iphop/`: host predictions from [iPHoP](https://bitbucket.org/srouxjgi/iphop) at genome level (`iphop_genome.csv`) and genus level (`iphop_genus.csv`), with Sankey plots of the predicted host taxonomy for each (full, and reduced to class level).
      - `spacepharer/`: CRISPR spacer-based host predictions from [SpacePHARER](https://github.com/soedinglab/spacepharer); only produced when `--predict_host_from_custom_spacers` is used.
        - `spacers_predictions.tsv`: matches between the provided spacers and the viral sequences.
        - `lineage_comparison.json`: rank-by-rank comparison between the predicted host and the source genome (MAG) the virus came from.
        - `custom_chosen_host.tsv`: host assigned to each virus from the provided spacer database.
        - `*sankey.html`: Sankey plots of the predicted host taxonomy.
    - `taxonomy/`
      - `genomad/`: [geNomad](https://github.com/apcamargo/genomad) taxonomy calls (`viruses.tsv`), taxonomy combined with metadata and counts (`viruses_genomad_taxonomy_counts.tsv`), interactive Krona and Sankey plots, and iTOL annotation files (`itol_genomad/`).
      - `viphogs/`: [ViPhOGs-based](https://www.mdpi.com/1999-4915/13/6/1164) taxonomy assignment (`viruses_reps_annotation.tsv`, `viruses_reps_annotation_taxonomy.tsv`), taxonomy tables (`viruses_modified*.tsv`, `viruses_viphogs_taxonomy_counts.tsv`), Krona and Sankey plots, and iTOL files (`itol_viphogs/`).
      - `vitap/`: [VITAP](https://github.com/DrKaiyangZheng/VITAP) taxonomy assignment and counts (`viruses_vitap_taxonomy_counts.tsv`), Krona and Sankey plots, and iTOL files (`itol_vitap/`).
    - `protein_clusters_phammseqs/`: proteins grouped into families ("phams") with [PhaMMseqs](https://github.com/chg60/phammseqs); only produced when `--phammseqs` is used.
      - `viruses_phams.tsv`: each protein (`proteinID`, with its description) and the pham it belongs to (`phamID`).

</details>

Viral sequences and prophages are pooled and clustered, and one representative per cluster is carried forward for annotation; its results are collected under `cluster_representatives/`. Taxonomy is assigned independently by three tools (geNomad, ViPhOGs and VITAP), each with its own directory under `taxonomy/` containing a lineage table, counts, and Krona, Sankey and iTOL visualisations. Hosts are predicted with iPHoP (genome and genus level) and, optionally, by CRISPR spacer matching with SpacePHARER. Lifestyle (virulent or temperate) is predicted with BACPHLIP. Proteins are annotated for AMR genes (AMRFinderPlus, DeepARG, RGI) and against the PVOG and TIGRFAM HMM databases, and these annotations are merged with the lifestyle call into `viruses_reps_annotated.gff.gz`.

## Plasmids

<details markdown="1">
<summary>Output files</summary>

- `plasmids/`
  - `plasmids-all-metadata.tsv.gz`: one row per **input** plasmid sequence (before clustering).
  - `plasmids.fasta.gz`, `plasmids.faa.gz`, `plasmids.gff.gz`: nucleotide sequences, proteins and GFF annotation of all plasmids (before clustering).
  - `plasmids_host_sankey.html`: Sankey plot of the taxonomy of the source genomes (MAGs) the plasmids were found in.
  - `cluster_representatives/`: results for the cluster representatives. Sequences are clustered with [vclust](https://github.com/refresh-bio/vclust) at 70% ANI and 50% coverage.
    - `plasmids.fasta.gz`, `plasmids.faa.gz`: nucleotide and protein sequences of one representative per plasmid cluster.
    - `plasmids.gff.gz`: GFF of the representatives, enriched with AMR and MOB typing results.
    - `functional_annotation/`
      - `mobility_stats_final.json`: summary counts combining plaSquid evidence and MOB-suite predicted mobility.
      - `amr/`
        - `plasmids.tsv`: antimicrobial resistance calls from [AMRFinderPlus](https://github.com/ncbi/amr), [DeepARG](https://github.com/gaarangoa/deeparg) and [RGI](https://github.com/arpcard/rgi).
        - `integrated_plasmids.gff`: GFF with the AMR annotations merged in.
      - `mobsuite_typer/`: results of [MOB-suite `mob_typer`](https://github.com/phac-nml/mob-suite#mob-typer).
        - `*_report.txt`: per-plasmid typing report, including predicted mobility.
        - `*biomarker_report.txt`: replicon, relaxase, mate-pair formation and oriT biomarkers found on each plasmid.
        - `*mge_report.txt`: other mobile genetic elements found on each plasmid.
      - `plasquid/`: predictions from [plaSquid](https://github.com/mgimenez720/plaSquid).
        - `final_report_per_contig.tsv`: plaSquid evidence summarised per plasmid.
        - `final_report_per_protein.tsv`: replication initiator protein (RIP) domain, MOB group and Inc group for each protein with a hit.
        - `mobility_classification.tsv`: mobility class assigned to each plasmid from plaSquid evidence.
        - `mobility_stats.json`: number of proteins with plaSquid hits, in total and for each evidence type.

</details>

Plasmid sequences are pooled and clustered, and one representative per cluster is carried forward for annotation; its results are collected under `cluster_representatives/`. The representatives are annotated for AMR genes (AMRFinderPlus, DeepARG, RGI), typed with MOB-suite and searched with plaSquid for replication initiator proteins, MOB groups and incompatibility groups.

## Catalogue summary

<details markdown="1">
<summary>Output files</summary>

- `catalogue.json`: summary counts for the whole run.

</details>

`catalogue.json` contains the following fields:

```text
{
  "total_sequences":              number of viral sequences, prophages and plasmids before QC filtering
  "unique_sequences":             number of viral sequences, prophages and plasmids after QC filtering
  "qc_excluded_sequences":        number of sequences removed at the QC step
  "viral_sequences":              number of viral sequences after QC
  "plasmids":                     number of plasmids after QC
  "prophages":                    number of prophages after QC
  "number_of_biomes":             number of unique biomes in the input samplesheet
  "viral_clusters":               number of viral clusters (= cluster representatives)
  "plasmid_clusters":             number of plasmid clusters (= cluster representatives)
  "total_proteins_viruses":       number of proteins predicted on viral sequences and prophages
  "total_proteins_plasmids":      number of proteins predicted on plasmids
  "total_proteins_viruses_reps":  number of proteins predicted on viral cluster representatives
  "total_proteins_plasmids_reps": number of proteins predicted on plasmid cluster representatives
}
```

## Additional data

Intermediate files produced while the inputs (from the samplesheet and, optionally, `--third_party_input`) are combined, quality-checked and split by sequence type. Files marked [`--publish_all`] are only published when the pipeline is run with `--publish_all`; all others are always published.

<details markdown="1">
<summary>Output files</summary>

- `additional_data/`
  - `rename_mapfile.tsv`: mapping between the original sequence names and the identifiers assigned by the pipeline.
  - `catalogue_filtered_metadata.tsv`: metadata for all sequences that passed QC filtering.
  - `catalogue_excluded_records_metadata.tsv`: records excluded for low quality, with the reason for exclusion.
  - `rename_contigs/`: input FASTA and GFF files combined across all samples, with unique catalogue-wide identifiers (`combined.fna`, `combined.gff`), and the mapping back to the original sequence names, biomes and types (`combined.tsv`).
  - `rename/` [`--publish_all`]: renamed `*.fna`, `*.faa` and `*.gff` files.
  - `viruses_checkv_report.tsv` [`--publish_all`]: [CheckV](https://bitbucket.org/berkeleylab/checkv) quality assessment for each sequence: completeness, contamination, provirus status and quality tier.
  - `barrnap/` [`--publish_all`]: rRNA gene predictions on the combined sequences (`combined_bac.gff`); skipped when `--skip_rrna_detection` is used.
  - `choose_sequences/` [`--publish_all`]: the sequences that passed QC filtering, with their metadata (`combined_filtered.tsv`, `combined_filtered.fna`, `combined_filtered.gff`).

</details>

Each input sequence is given a stable catalogue identifier (`MGYV*` by default; configurable with `--start_accession`/`--end_accession`), so that downstream steps and file names do not depend on the original sample naming; `rename_mapfile.tsv` maps these back to the original names. CheckV and barrnap are run on the renamed sequences. The quality filter then removes low-quality sequences (for example, sequences with no detected viral genes), listed in `catalogue_excluded_records_metadata.tsv`, and splits the remaining sequences into the virus/prophage and plasmid sets used by the [Viruses](#viruses) and [Plasmids](#plasmids) steps above.

## MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML report that can be viewed in a web browser.
  - `multiqc_data/`: parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) generates a single HTML report summarising the whole run. Most of the pipeline QC results are shown in the report, and further statistics are available in the report data directory.

The report collates QC results from supported tools such as CheckV, and also lists the software versions used, for traceability. For more information about MultiQC reports, see <http://multiqc.info>.

## Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reformatted samplesheet used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used for the run: `params.json`.
  - Collated software versions: `metaviraverse_software_mqc_versions.yml`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) generates reports on how the pipeline ran. Use them to troubleshoot failures and to check launch commands, run times and resource usage.
