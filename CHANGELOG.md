# EBI-Metagenomics/viraverse: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [date]

Initial release of EBI-Metagenomics/metaviraverse, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- `catalogue_manifest.json`: accession range, deduplication and clustering settings, software versions and the checksummed list of files that a catalogue update reads (first step towards `--mode update`).
- `--mode generate|update` and `--existing_catalogue` (in development): in update mode the existing catalogue is read through its manifest and every listed file is checked (size and SHA-256) before the run continues; the report is published as `update/existing_catalogue_check.json`. New sequences are not yet merged into the existing catalogue.
- With `--viral_clustering vjdb`, the vjdb representative sequences are published (`viruses/clustering/viruses_vjdb_reps.fna.gz`), the membership table is renamed to `viruses/clustering/viruses_vjdb_clusters.csv`, and the run stops if vjdb-clustering finds identical sequences that `CHOOSE_SEQUENCES` kept.

### `Fixed`

### `Dependencies`

### `Deprecated`
