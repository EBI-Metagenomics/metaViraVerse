# EBI-Metagenomics/viraverse: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [date]

Initial release of EBI-Metagenomics/metaviraverse, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- `--analyses_metadata` (optional): metadata of MGnify analyses (`analyses_metadata.tsv` from `scripts/collect_mobilome_from_analyses.py`: accession, study_id, sample_id, biome_id). Sequences named `<ERZ accession>|...` get `Source_accession`, `Source_sample` and `Source_project` (and `Biome` when missing) in `viruses-all-metadata.tsv.gz` and `plasmids-all-metadata.tsv.gz`.
- `README.md` of the catalogue in `--outdir` (`GENERATE_README`, `bin/generate_catalogue_readme.py`): summary numbers, tools and annotations that were run with their versions, how the catalogue was built, metadata statistics (counts by `Source_of_prediction`, `Source` and `checkv_quality`, completeness and contamination ranges, missing values per column) and links to the documentation. New optional params for its title section: `--catalogue_name`, `--catalogue_version`, `--catalogue_website_url`, `--source_catalogues`.

- Tables behind the Sankey plots (`*_sankey.tsv`, `--output-table` of `bin/plot_taxonomy_sankey.py`) are published next to each plot: taxonomy (`genomad_`, `viphogs_`, `vitap_sankey.tsv`), hosts (`iphop_genome_`, `iphop_genus_`, `spacepharer_` with their `_reduced` versions) and source genomes (`viruses_source_`, `plasmids_source_sankey.tsv`). One row per link: source and target node, their levels, the target lineage and the number of sequences.

### `Fixed`

- `catalogue.json` sequence counts: `total_sequences` is now the number of input sequences (rows of the rename map; it was the number after deduplication), and `qc_excluded_sequences` no longer includes duplicates. New keys `duplicate_sequences`, `deduplicated_sequences`, `rrna_excluded_sequences` and `not_determined_excluded_sequences` (also in `choose_sequences.py`'s stats JSON as `number_of_unique_sequences` and `number_of_<reason>_excluded`). The README summary and `scripts/check_catalogue_results.py` use them; the check script also verifies that the counts add up.
- `collect_catalogue_stats.py`: `viral_clusters` and `plasmid_clusters` in `catalogue.json` were one lower than the real number. One line was subtracted for a header, but the reps list it reads has none. A header is now skipped only when present.
- `collect_metadata.py`: sequences of MGnify analyses (`ERZ..._N|type-start:end`, no space in the name) were not matched to their catalogue ID, so `ID`, `Cluster_rep`, taxonomy and hosts were missing and `Original_ID` had `|` replaced by `-`. Rows are now joined on the catalogue ID. `--analyses_metadata` is also joined on the ERZ accession without the contig number.

### `Dependencies`

### `Deprecated`
