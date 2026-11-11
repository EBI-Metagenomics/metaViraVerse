# EBI-Metagenomics/viraverse: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [date]

Initial release of EBI-Metagenomics/metaviraverse, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- `--analyses_metadata` (optional): metadata of MGnify analyses (`analyses_metadata.tsv` from `scripts/collect_mobilome_from_analyses.py`: accession, study_id, sample_id, biome_id). Sequences named `<ERZ accession>|...` get `Source_accession`, `Source_sample` and `Source_project` (and `Biome` when missing) in `viruses-all-metadata.tsv.gz` and `plasmids-all-metadata.tsv.gz`.

### `Fixed`

### `Dependencies`

### `Deprecated`
