# EBI-Metagenomics/viraverse: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [date]

Initial release of EBI-Metagenomics/metaviraverse, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- `--analyses_metadata` (optional): metadata of MGnify analyses (`analyses_metadata.tsv` from `scripts/collect_mobilome_from_analyses.py`: accession, study_id, sample_id, biome_id). Sequences named `<ERZ accession>|...` get `Source_accession`, `Source_sample` and `Source_project` (and `Biome` when missing) in `viruses-all-metadata.tsv.gz` and `plasmids-all-metadata.tsv.gz`.
- `README.md` of the catalogue in `--outdir` (`GENERATE_README`, `bin/generate_catalogue_readme.py`): summary numbers, tools and annotations that were run with their versions, how the catalogue was built and links to the documentation. New optional params for its title section: `--catalogue_name`, `--catalogue_version`, `--catalogue_website_url`, `--source_catalogues`.

### `Fixed`

- `collect_metadata.py`: sequences of MGnify analyses (`ERZ..._N|type-start:end`, no space in the name) were not matched to their catalogue ID, so `ID`, `Cluster_rep`, taxonomy and hosts were missing and `Original_ID` had `|` replaced by `-`. Rows are now joined on the catalogue ID. `--analyses_metadata` is also joined on the ERZ accession without the contig number.

### `Dependencies`

### `Deprecated`
