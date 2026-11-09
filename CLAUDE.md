# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

metaViraVerse is an nf-core-style Nextflow DSL2 pipeline (EBI MGnify) that builds a **viral catalogue** and a **plasmid catalogue** from already-predicted viral sequences, prophages and plasmids. Input is MGnify MAG catalogue results (VIRify + mobilome-annotation-pipeline output, `--input`) and/or third-party FASTA/GFF (`--third_party_input`). User docs: `docs/mgnify_usage.md`, `docs/third_party_usage.md`, `docs/methods.md` (method details), `docs/output.md` (output layout and metadata column definitions).

## Commands

```bash
# Python unit tests for bin/ scripts (CI: .github/workflows/pytest.yml, Python 3.11–3.13)
pip install -r requirements-test.txt pytest
pytest tests/unit                                   # all
pytest tests/unit/test_rename_contigs.py            # one file
pytest tests/unit/test_rename_contigs.py -k TestDefinePrefix   # one test class/case

# nf-test (module/subworkflow/pipeline tests; config in nf-test.config, profile "test")
nf-test test modules/local/sankey_plot/tests/main.nf.test --profile +docker
nf-test test tests/default.nf.test --profile +docker

# Lint (CI: .github/workflows/linting.yml)
pre-commit run --all-files                          # prettier, trailing whitespace, EOF
nf-core pipelines lint

# Run
nextflow run main.nf -profile <profile> --outdir <dir> --input samplesheet.csv --catalogues_metadata genomes-all_metadata.tsv
```

Unit tests import scripts directly from `bin/` (`sys.path.insert(0, BIN_DIR)`); fixtures live in `tests/unit/fixtures/<script>/`. `.gitignore` has a blanket `*.tsv` rule, so new `.tsv` test fixtures must be force-added (`git add -f`).

## Code layout

- `workflows/metaviraverse.nf` — top level: `THIRD_PARTY_DATA` (optional) → `PREPROCESSING` → `PROCESS_VIRAL_SEQUENCES` and `PROCESS_PLASMIDS` → `COLLECT_CATALOGUE_STATS` (`catalogue.json`) and `COLLECT_METADATA` (`*-all-metadata.tsv.gz`, `viruses-cluster-stats.tsv.gz`).
- `subworkflows/local/` — pipeline logic; `modules/local/` — thin wrappers that mostly call a script in `bin/`; `modules/nf-core/` — vendored nf-core modules (don't edit; tracked in `modules.json`).
- `bin/` — Python (and plaSquid R) scripts run inside processes. `scripts/` — helper scripts run **outside** the pipeline on EBI NFS to prepare inputs (`collect_data_from_catalogues.py` builds the samplesheet; `collect_crispr_spacers_from_catalogues.py` builds `--custom_spacers_fasta/metadata`).
- `conf/modules.config` — all `ext.args`/`ext.prefix`. Publishing uses workflow outputs: the `output` block in `main.nf` defines every published path and file name (no `publishDir`); `outputDir = params.outdir` in `nextflow.config`.
- `development/` — local scratch data and HTML report prototypes; not part of the pipeline.

## Data flow and identifiers (needs several files to see)

1. **Rename** (`bin/rename_contigs.py`): every input sequence gets a catalogue ID `<rename_accession><n>` (`MGYV` + 10-digit zero-padded for MGnify). It writes the rename map with columns `original, temporary, short, biome, type, definition`: `temporary` = new ID, `short` = original ID up to the first space, `type` = `genome`/`metagenome`, `definition` = `virus`/`prophage`/`plasmid` (prefixed `third_party_` for third-party rows). Downstream code maps IDs back through this file.
2. **Separate** (`bin/separate_sequences.py`) into virus/prophage/plasmid using the map's `definition`. CheckV and barrnap run on viruses and prophages only.
3. **Choose/deduplicate** (`bin/choose_sequences.py`): deduplicates across all three categories by sequence SHA-256, with priority MGnify-metagenome > MGnify-genome > third-party-metagenome > third-party-genome; biomes of duplicates are merged. Low-quality records go to the excluded table.
4. **Viruses** = viral sequences + prophages pooled, clustered with vclust (95% ANI / 85% coverage; `conf/modules.config` `VCLUST_*`) or, with `--viral_clustering vjdb`, with the included pipeline `pipelines/vjdb-clustering` (git submodule; its membership CSV is converted to the vclust `object`/`cluster` TSV by `VJDB_CLUSTERS_TSV`), then **only cluster representatives** are annotated: taxonomy (geNomad, VITAP, ViPhOGs), host (iPHoP; SpacePHARER with custom spacers), lifestyle (BACPHLIP), proteins (HMMER, AMR, PhaMMseqs). `BUILD_FINAL_GFF` merges annotations into one GFF. Protein steps live in `subworkflows/local/proteins_subwf` (`PROTEINS_PROCESSING`), which is called from both the viral and the plasmid subworkflow with different skip flags.
5. **Plasmids** clustered separately (vclust: pairs with ≥70% ANI / ≥50% coverage, Leiden at gANI 0.35), representatives annotated with AMR, MMseqs2 protein clustering, MOB-suite `mob_typer` and plaSquid (R scripts); `bin/annotate_plasmid_gff.py` merges these into the GFF; `bin/update_mobility_stats.py` combines plaSquid and MOB-suite mobility counts.

Original MGnify sequence names look like `MGYG000519004_12|viral_sequence-1:32972` (genome, contig, type, coordinates); the source genome is the part before the first `_`, which is how scripts join to `--catalogues_metadata` (lineage, sample, study).

## File formats

Identifiers: after renaming, **sequences** are catalogue IDs (`MGYV0000000001`), but **proteins keep their source IDs** (`MGYG000999448_03173`, with a description in FASTA headers). A protein's representative cannot be derived from its ID; use the reps GFF (CDS rows lie on the representative's seqid). All files below may be plain or gzip/bgzip; bgzipped published files (`*.gz` with `.tbi`/`.fai`/`.gzi`) read with plain gzip.

Sequences and annotation

- Renamed GFF / reps GFF (`viruses_reps.gff.gz`, published sorted, bgzipped, tabix-indexed): GFF3, `##gff-version 3`; seqid (col 1) = catalogue ID; one region feature per sequence (source `VIRify`, type `viral_sequence`/`prophage`/`plasmid`, `ID=<MGYV…>` plus VIRify/CheckV attributes), then `CDS` rows with `ID=<source protein ID>`. No `##FASTA` section.
- `viruses_reps_annotated.gff.gz` (`BUILD_FINAL_GFF`): the reps GFF with BACPHLIP (`bacphlip_virulent`/`bacphlip_temperate` on the region), HMMER (`<db>_accession`, … on CDS) and AMR attributes added.
- `integrated_viruses.gff` (AMRINTEGRATOR): GFF of AMR hits, seqid = catalogue ID.
- Reps FASTA (`viruses_reps.fasta.gz` / `.faa.gz`): FNA headers = catalogue ID; FAA headers = source protein ID + description.

Clustering

- vclust cluster TSV (internal `*_clusters.tsv`): header `object\tcluster`, one row per sequence, `cluster` = representative ID. The blastn/CheckV variant has no header (`rep\tmembers`).
- vjdb membership CSV (`viruses_vjdb_clusters.csv`): `vjdb_id,copy_of,seq_gc,seq_len,seq_hash,cluster_rep`; `seq_hash` = case-sensitive MD5 of the sequence; `copy_of` must equal `vjdb_id` (identical sequences were removed before). `viruses_vjdb_reps.fna.gz`: one sequence per `cluster_rep`.
- `--mode update` only: `viruses_new_sequences_clusters.csv` (vjdb `*_new_clusters.csv`: the CSV columns plus `status` = `new_unique`/`new_copy`/`existing_copy`/`resubmitted`/`duplicate_in_batch`) and `viruses_cluster_events.tsv` (`cluster_rep, event, old_cluster_reps, n_existing, n_new`; `event` = `new`/`grown`/`merged`; `old_cluster_reps` `;`-separated; unchanged clusters are not listed). Internal `*_clusters_to_annotate.tsv` (vclust format, `new` clusters only) and `*_dropped_reps.txt` (one old representative per line, absorbed by a merge) come from `bin/vjdb_clusters_tsv.py`.

Per-representative annotation tables (ID column → what it holds)

- geNomad `viruses.tsv`: TSV, `seq_name` (catalogue ID), …, `taxonomy` (`;`-separated).
- VITAP best: TSV, `Genome_ID` (catalogue ID), `lineage`.
- ViPhOGs taxonomy `viruses_reps_annotation_taxonomy.tsv`: TSV, header `contig_ID` (with leading spaces: strip header names), then rank columns `superkingdom … genus`; all-empty ranks = no assignment.
- ViPhOGs annotation `viruses_reps_annotation.tsv`: TSV `Contig` (catalogue ID), `CDS_ID` (source protein ID), `Start, End, Direction, Best_hit, Abs_Evalue_exp, Label`. `viruses_modified*.tsv`: `ViPhOG, query` (protein ID), `Ratio, Abs_Evalue_exp, Taxon`.
- iPHoP `iphop_genome.csv` / `iphop_genus.csv`: CSV, `Virus` (catalogue ID), host lineage in `Host taxonomy` / `Host genus`, `Confidence score`; several rows per virus (best = highest score); chunk concatenation repeats the header.
- SpacePHARER `custom_chosen_host.tsv`: TSV, `viral_seq` (catalogue ID), `host_lineage` (`NA` = none), `evalue`; several rows per virus (best = lowest evalue).
- BACPHLIP `*.bacphlip`: TSV, header `\tVirulent\tTemperate` (first column unnamed = catalogue ID); published bgzipped + tabix (`viruses.bacphlip.gz`).
- HMMER `*.tbl.gz`: `hmmsearch --tblout`, whitespace-separated, `#` comment lines; col 1 = target = source protein ID, col 3 = profile accession. Several hits per protein (`BUILD_FINAL_GFF` keeps the first). One file per HMM database, named `<meta.id>_<db>.tbl.gz`; `*_summary.tsv.gz` per database.
- AMR: `viruses_amrfinderplus.tsv` (AMRFinderPlus 4 TSV, protein ID in `Protein id`); `viruses_deeparg.tsv`, `viruses_rgi.tsv` (hAMRonization TSV, protein ID in `input_sequence_id`).
- PhaMMseqs `viruses_phams.tsv`: `phamID\tproteinID function` (protein ID + description); whole-catalogue analysis, not per representative.

Metadata tables

- Combined/filtered/catalogue metadata (`CHOOSE_SEQUENCES`; `catalogue_filtered_metadata.tsv`): `sequence_id, original_name, description, type, source_of_prediction, biomes, sequence_length, rrna, sequence_sha256`, the CheckV `QUALITY_COLUMNS` (`provirus … warnings`), `definition`. `sequence_sha256` = SHA-256 of the upper-cased, stripped sequence. The excluded table adds `filter_reason` (quality filters, `duplicate`, `in_existing_catalogue`).
- `viruses-all-metadata.tsv.gz` (`COLLECT_METADATA`): one row per input viral sequence: `ID, Original_ID, Source_of_prediction, rRNA, Cluster_rep, vitap_lineage, viphogs_lineage, genomad_lineage, Host_iphop_genome_lineage, Host_iphop_genus_lineage, Host_spacepharer_lineage, Source, Biome, Source_accession, Source_lineage, Source_sample, Source_project, Sequence_length, Sequence_sha256` + CheckV fields.
- `viruses-cluster-stats.tsv.gz`: one row per representative: `Rep_ID, Rep_vitap_lineage, Rep_viphogs_lineage, Rep_genomad_lineage, Rep_source, Rep_biome, cluster_size, cluster_mean_viral_genes, cluster_mean_gene_count, cluster_biomes, cluster_types, cluster_completeness_range, cluster_contamination_range`.
- `--catalogues_metadata` (MGnify `genomes-all_metadata.tsv`): joined on `Genome` (the part of the original name before the first `_`); uses `Genome_accession, Lineage, Sample_accession, Study_accession`.

Catalogue description and update

- `catalogue_manifest.json` (`bin/write_catalogue_manifest.py`): `manifest_version` (1), `catalogue` (`created, mode, parent, pipeline_version, run_name, commit`), `accessions` (`prefix, first, last, count, zero_padded_digits`), `deduplication`, `viral_clustering` (`method, update_ready, not_update_ready_reason`), `software_versions`, `files`: key → `{path, size, sha256}` (a list for `viruses.hmmer_tables`); `path` is relative to `--outdir`. Keys are written by `manifestInput()` in `workflows/metaviraverse.nf` and read by `loadExistingCatalogue()` in `subworkflows/local/existing_catalogue/main.nf`: keep them in sync.
- `update/existing_catalogue_check.json` (`bin/check_existing_catalogue.py`): `status, manifest_sha256, manifest_version, catalogue, accessions, files_checked`.

## Conventions

- **Optional inputs/steps** (`--skip_*`, `--predict_host_from_custom_spacers`, `--annotation_db`): skipped steps emit `channel.empty()`. Before passing such a channel into a process, map to the file and add `.ifEmpty([])`; in the module build the flag conditionally (`def x_arg = x ? "--x ${x}" : ""`); in the script make the argument `required=False` and handle its absence. Otherwise the downstream process never runs.
- After `.join(..., remainder: true)` the `.map { ... }` closure must take exactly as many values as the joined tuple has; update it whenever a join is added.
- `FIND_CONCATENATE` (used to merge chunked outputs) names its output exactly `ext.prefix`, with no extension unless the prefix has one — publish targets that keep the source name (ending in `/`) inherit that.
- To publish a new file, carry it up in the subworkflow's emitted record (`Preprocessed`, `ViralCatalogue`, `PlasmidCatalogue`, `CatalogueSummary`, and the nested `TaxonomyResult`, `HostResult`, `ProteinsResult`, …) and add a `field >> "target"` line to the `output` block. A null field is skipped, but the target expression is still evaluated, so use `?.` on nested optional records. Two files sent to the same target fail the run.
- Intermediate files are published only with `--save_intermediates` (the `*_intermediates` outputs, `enabled params.save_intermediates`); keep `docs/output.md` in sync with which outputs are always published.
- Missing values in final metadata tables are written as the string `missing` (upstream steps may write `NA` / `not-provided`).
- Scripts must read both plain and gzip-compressed inputs; many intermediate files are `.gz`.
- Chunked steps (CheckV, iPHoP, MOB-suite, VITAP) concatenate chunk outputs, which can repeat header lines — parsers should skip repeated headers.
