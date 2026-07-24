# EBI-Metagenomics/metaviraverse: Output

## Introduction

This document describes the output produced by the pipeline. The pipeline builds a viral catalogue and a plasmid catalogue from predicted viral sequences, prophages, and plasmids, and organises the results under two top-level groups, **viruses/** and **plasmids/**, plus supporting directories for intermediate data, QC, and pipeline metadata.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory. An example output tree (from a run on barley rhizosphere samples) is used throughout to illustrate real file names.

## Pipeline overview

The pipeline processes data using the following steps:

- [Catalogue summary](#catalogue-summary) - top-level statistics and combined metadata for the whole catalogue
- [Viruses](#viruses) - clustered viral sequences and prophages, with taxonomy, host, lifestyle and functional annotation
- [Plasmids](#plasmids) - clustered plasmid sequences
- [Additional data](#additional-data) - intermediate files produced while preparing the input (renaming, quality control, rRNA detection, sequence splitting)
- [MultiQC](#multiqc) - aggregate report describing results and QC from the whole pipeline
- [Pipeline information](#pipeline-information) - report metrics generated during the workflow execution

## Catalogue summary

<details markdown="1">
<summary>Output files</summary>

- `catalogue.json`: summary counts for the whole run (input/unique/excluded sequences, viral sequences, prophages, plasmids, number of biomes, viral/plasmid cluster counts, and total protein counts for viruses and plasmids).
- `viruses-all-metadata.tsv.gz`: one row per **input** viral sequence/prophage (before clustering), combining source, quality, and taxonomy information.
- `plasmids-all-metadata.tsv.gz`: one row per **input** plasmid sequence (before clustering).
- `viruses-cluster-stats.tsv.gz`: one row per viral **cluster representative**, with per-cluster statistics (cluster size, mean gene/viral-gene counts, biomes and source types in the cluster, completeness/contamination ranges).

</details>

`catalogue.json` gives an at-a-glance summary of the run, e.g.:

```json
{
  "total_sequences": 1283,
  "unique_sequences": 1253,
  "qc_excluded_sequences": 30,
  "viral_sequences": 141,
  "plasmids": 1057,
  "prophages": 55,
  "number_of_biomes": 1,
  "viral_clusters": 195,
  "plasmid_clusters": 1008,
  "total_proteins_viruses": 15191,
  "total_proteins_plasmids": 25316
}
```

`viruses-all-metadata.tsv.gz` and `plasmids-all-metadata.tsv.gz` list every sequence that went into the catalogue (identifiers, originating sample/project/biome, taxonomic lineage of the source genome, sequence length/checksum, and — for viruses — CheckV quality metrics and the taxonomy assigned by VITAP, ViPhOGs and geNomad). `viruses-cluster-stats.tsv.gz` summarises each viral cluster by its representative sequence.

## Viruses

Viral sequences and prophages predicted upstream (e.g. by [VIRify](https://github.com/EBI-Metagenomics/emg-viral-pipeline) or [mobilome-annotation-pipeline](https://github.com/EBI-Metagenomics/mobilome-annotation-pipeline)) are pooled into a single **viruses** set, clustered at 95% ANI, and the cluster representatives are annotated.

<details markdown="1">
<summary>Output files</summary>

- `viruses/`
  - `viruses.fasta.gz`: all viral sequences and prophages, combined, before clustering.
  - `cluster_reps/`
    - `viruses.fasta.gz`, `viruses.faa.gz`, `viruses.gff.gz` (+ `.fai`/`.gzi`/`.tbi` indices): nucleotide, protein, and GFF records for one representative sequence per viral cluster.
    - `viruses_final.gff.gz`: the representative GFF enriched with lifestyle (BACPHLIP), HMMER (e.g. PVOG, TIGRFAM) and AMR annotations in a single file.
    - `functional_annotation/`
      - `amr/`: antimicrobial resistance calls from AMRFinderPlus, DeepARG and RGI (`viruses.tsv`), plus the GFF with AMR annotations merged in (`integrated_viruses.gff`).
      - `hmmer/`: HMMER hits against hmm profile databases (`*.tbl.gz`) with per-accession function summaries (`*_summary.tsv.gz`).
    - `host_detection/`
      - `iphop/`: host predictions from [iPHoP](https://bitbucket.org/srouxjgi/iphop) at genome level (`iphop_genome.csv`) and genus level (`iphop_genus.csv`).
      - `spacepharer/`: CRISPR spacer-based host predictions from SpacePHARER (`spacers_predictions.tsv`); only produced when `--predict_host_from_custom_spacers` is used.
    - `lifestyle/`: virulent/temperate lifestyle prediction from BACPHLIP (`viruses.bacphlip.gz`).
    - `taxonomy/`
      - `genomad/`: [geNomad](https://github.com/apcamargo/genomad) taxonomy calls (`viruses.tsv`), taxonomy combined with metadata and counts (`viruses_genomad_taxonomy_counts.tsv`), interactive Krona/Sankey plots, and iTOL annotation files (`itol_genomad/`).
      - `viphogs/`: ViPhOGs-based taxonomy assignment (`viruses_reps_annotation.tsv`, `viruses_reps_annotation_taxonomy.tsv`), taxonomy tables (`viruses_modified*.tsv`, `viruses_viphogs_taxonomy_counts.tsv`), Krona/Sankey plots, and iTOL files (`itol_viphogs/`).
      - `vitap/`: [VITAP](https://github.com/DrKaiyangZheng/VITAP) taxonomy assignment and counts (`viruses_vitap_taxonomy_counts.tsv`), Krona/Sankey plots, and iTOL files (`itol_vitap/`).
    - `phammseqs/`: protein sequences grouped into phamilies by MMseqs2 clustering; only produced when `--phammseqs` is enabled.

</details>

Viral sequences and prophages are first pooled (`viruses.fasta.gz`) and clustered at 95% average nucleotide identity; one representative per cluster is carried forward for downstream annotation, and its results are collected under `cluster_reps/`. Taxonomy is assigned independently by three tools (geNomad, ViPhOGs, VITAP), each contributing its own directory under `taxonomy/` with a lineage table, counts, and Krona/Sankey/iTOL visualisations. Host predictions come from iPHoP (genome- and genus-level) and, optionally, from CRISPR spacer matching via SpacePHARER. Lifestyle (virulent vs. temperate) is predicted with BACPHLIP. Predicted proteins are functionally annotated for AMR determinants (AMRFinderPlus, DeepARG, RGI) and against the PVOG/TIGRFAM HMM databases, and folded back into `viruses_final.gff.gz` alongside the lifestyle call.

## Plasmids

<details markdown="1">
<summary>Output files</summary>

- `plasmids/`
  - `plasmids.fasta.gz`: all input plasmid sequences, combined, before clustering.
  - `cluster_reps/`
    - `plasmids.fasta.gz`, `plasmids.faa.gz`, `plasmids.gff.gz` (+ `.fai`/`.gzi`/`.tbi` indices): nucleotide, protein, and GFF records for one representative sequence per plasmid cluster.

</details>

Plasmid sequences are pooled (`plasmids.fasta.gz`) and clustered at 35% (default; configurable) global ANI. Only cluster-representative sequences, proteins, and coordinates are retained under `cluster_reps/` — plasmids do not currently go through taxonomy, host detection, or functional annotation (host detection and mobility-element typing are planned).

## Additional data

Intermediate files produced while the raw inputs (from the samplesheet, and optionally `--third_party_input`) are combined, quality-checked, and split by sequence type.

<details markdown="1">
<summary>Output files</summary>

- `additional_data/`
  - `rename_contigs/`: input FASTA/GFF combined across all samples with unique catalogue-wide identifiers (`combined.fna`, `combined.gff`), and the mapping back to original sequence names, biomes and types (`combined.tsv`).
  - `quality/`: [CheckV](https://bitbucket.org/berkeleylab/checkv) quality assessment per sequence (`combined.tsv`, `quality_summary.tsv`) — completeness, contamination, provirus status, and quality tier.
  - `barrnap/`: rRNA gene predictions on the combined sequences (`combined_bac.gff`); skipped when `--skip_rrna_detection` is set.
  - `choose_sequences/`: sequences split by type after quality filtering — combined metadata for all input sequences (`combined_metadata.tsv`), the metadata that passed filtering (`combined_filtered.tsv`, `combined_filtered.fna`, `combined_filtered.gff`), and the records excluded for low quality with their exclusion reason (`combined_excluded.tsv`).

</details>

`rename_contigs/` assigns each incoming sequence a stable catalogue identifier (`MGYV*` by default, configurable via `--start_accession`/`--end_accession`) so that downstream steps and file names are independent of the original sample naming. `quality/` and `barrnap/` run on this renamed set. `choose_sequences/` then applies the quality filter (excluding, for example, sequences with no detected viral genes) and splits the remaining sequences into the viral/prophage and plasmid pools that feed the [Viruses](#viruses) and [Plasmids](#plasmids) workflows above.

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. CheckV. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.
  - Collated software versions: `metaviraverse_software_mqc_versions.yml`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
