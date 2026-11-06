# EBI-Metagenomics/viraverse: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [date]

Initial release of EBI-Metagenomics/metaviraverse, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- `catalogue_manifest.json`: accession range, deduplication and clustering settings, software versions and the checksummed list of files that a catalogue update reads (first step towards `--mode update`).
- `--mode generate|update` and `--existing_catalogue` (in development): in update mode the existing catalogue is read through its manifest and every listed file is checked (size and SHA-256) before the run continues; the report is published as `update/existing_catalogue_check.json`. In update mode new accessions continue after the existing catalogue, sequences the existing catalogue already has keep their ID (excluded as `in_existing_catalogue`, biomes added to the existing record), and the rename map and catalogue metadata tables cover the whole catalogue. In update mode new viruses are clustered with the existing catalogue's representatives (vjdb-clustering update; stable cluster IDs, cluster events published in `viruses/clustering/update/`) and only the representatives of new clusters are annotated; annotations are not yet merged with the existing ones, so an updated catalogue is marked not update-ready.
- With `--viral_clustering vjdb`, the vjdb representative sequences are published (`viruses/clustering/viruses_vjdb_reps.fna.gz`), the membership table is renamed to `viruses/clustering/viruses_vjdb_clusters.csv`, and the run stops if vjdb-clustering finds identical sequences that `CHOOSE_SEQUENCES` kept.

### `Fixed`

- The RGI database is prepared once per run (`PREPARE_RGI_DB`) and shared by the viral and the plasmid AMR annotation; before, CARD was downloaded and processed separately for viruses and for plasmids. `--rgi_db` may now also be an unprocessed directory with `card.json`, which is processed once.
- Stub runs without containers: the RGI module stubs no longer call `rgi`, `seqkit/split2` has a stub, and the tabix version is read with `sed` instead of `grep -P` (not available in BSD grep).

### `Dependencies`

### `Deprecated`
