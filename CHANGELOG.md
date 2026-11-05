# EBI-Metagenomics/viraverse: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [date]

Initial release of EBI-Metagenomics/metaviraverse, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- `catalogue_manifest.json`: accession range, deduplication and clustering settings, software versions and the checksummed list of files that a catalogue update reads (first step towards `--mode update`).
- With `--viral_clustering vjdb`, the vjdb representative sequences are published (`viruses/clustering/viruses_vjdb_reps.fna.gz`), the membership table is renamed to `viruses/clustering/viruses_vjdb_clusters.csv`, and the run stops if vjdb-clustering finds identical sequences that `CHOOSE_SEQUENCES` kept.

### `Fixed`

### `Dependencies`

### `Deprecated`
