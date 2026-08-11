# Third-party data usage

## Overview

Third-party input lets sequences from any external source — a public database dump, a collaborator's dataset, an in-house collection, a previously published catalogue — be folded into the same metaViraVerse run as MGnify-derived data, so that both end up in one combined, deduplicated catalogue. Unlike MGnify input, third-party records don't need to already carry GFF/protein annotations: only a nucleotide FASTA is required, and the pipeline calls genes for you.

Third-party input is independent of MGnify input: a run can supply MGnify data only, third-party data only, or both together. See [MGnify usage](mgnify_usage.md) for the MGnify-side samplesheet, and [Methods](methods.md) for the full, step-by-step description of what happens to every sequence once it's in the pipeline.

## Input samplesheet

Third-party data is provided via `--third_party_input`, a CSV samplesheet with one row per FASTA file, following [`assets/schema_third_party_input.json`](../assets/schema_third_party_input.json):

```csv
id,fasta,type,source,biome
my_viruses,/path/to/my_viruses.fasta.gz,virus,metagenome,human gut
my_plasmids,/path/to/my_plasmids.fasta,plasmid,genome,marine sediment
```

| Column             | Required | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
| ------------------ | -------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `id`               | Yes      | Unique identifier for the record.                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
| `fasta`            | Yes      | FASTA file with the nucleotide sequence(s) for the chosen `type`. May be gzip-compressed.                                                                                                                                                                                                                                                                                                                                                                                               |
| `type`             | Yes      | `virus` (free/lytic viral sequences), `prophage`, or `plasmid` — the sequence's biological category.                                                                                                                                                                                                                                                                                                                                                                                    |
| `source`           | No       | `genome` (the sequence comes from a MAG or isolate) or `metagenome` (the sequence comes from an assembly) — the same distinction, and the same column name, as the MGnify samplesheet's own `source` column. Used later to decide which copy of a duplicated sequence to keep; if omitted, the sequence is still processed, but it is given the _lowest_ deduplication priority of all (lower than an explicit `genome`), so we recommend setting it whenever the true origin is known. |
| `biome`            | No       | Free-text environmental biome, ideally no more than two words (e.g. `marine sediment`, `human gut`).                                                                                                                                                                                                                                                                                                                                                                                    |
| `study_accession`  | No       | INSDC study accession (primary or secondary) associated with the data, kept as metadata.                                                                                                                                                                                                                                                                                                                                                                                                |
| `sample_accession` | No       | INSDC sample accession associated with the data, kept as metadata.                                                                                                                                                                                                                                                                                                                                                                                                                      |

A minimal samplesheet needs only `id`, `fasta` and `type` — everything else is optional metadata that is carried through to the final catalogue but doesn't change how sequences are processed (with the exception of `source`, see [Integration with the MGnify catalogue](#integration-with-the-mgnify-catalogue) below).

## Preparation steps

Before third-party sequences are merged with anything else, each input FASTA goes through its own short preparation subworkflow, since — unlike MGnify input — it can't be assumed to already be well-formed or gene-called:

1. **Decompression** — gzip-compressed FASTA files are detected by their `.gz` extension and decompressed; uncompressed files pass through unchanged.
2. **Validation** — every FASTA is checked with [FALint](https://github.com/GallVp/fa-lint), a FASTA linter/validator. Records from a file that fails validation are **not** carried forward into the rest of the pipeline (they are dropped rather than causing the whole run to fail); every invalid file is logged, together with its sample ID, to `invalid_fastas.tsv` in the output directory, so that malformed input can be noticed and fixed without having to re-read pipeline logs.
3. **Gene calling** — validated sequences are run through [Pyrodigal](https://github.com/althonos/pyrodigal) to predict genes, producing a GFF of CDS features and the corresponding translated protein FASTA. This is what lets a third-party record — which may have arrived as nothing more than raw nucleotide sequence — carry the same `fasta`/`gff`/`faa` triple that MGnify records already have, so that from this point onward the two sources are indistinguishable in shape and can be processed by exactly the same downstream steps.

## Integration with the MGnify catalogue

Once prepared, third-party records are merged with any MGnify input at the very first step of pre-processing (renaming) and flow through the same pipeline from then on — see [Methods: Pre-processing](methods.md#pre-processing) for the full description. The parts specific to how third-party data is integrated are:

- **Accession numbering** — MGnify sequences are renamed first; third-party sequences are renamed immediately afterwards, continuing the same accession counter rather than starting a new range. The combined catalogue therefore has one contiguous numbering range regardless of how many records came from each source.
- **Provenance tracking** — the rename step's mapping file records, for every sequence, a `type` column (`genome`/`metagenome`, taken from MGnify's own `source` column or the third-party samplesheet's `source` column) and a `definition` column (the sequence's biological category). For third-party records, `definition` is copied from the samplesheet's `type` column but prefixed with `third_party_` (e.g. `third_party_prophage`) — that prefix is what lets every later step tell an MGnify sequence and a third-party sequence apart, even once both are described using the same `genome`/`metagenome`/`virus`/`prophage`/`plasmid` vocabulary.
- **Separation** — sequences are split into virus/prophage/plasmid groups using that `definition` column, with the `third_party_` prefix stripped for the purpose of bucketing: a `third_party_prophage` sequence is filed into the same group as a plain `prophage` one, so MGnify and third-party sequences of the same category are processed together from here on (clustering, quality evaluation, etc. all operate on the merged group, not on MGnify and third-party separately).
- **Deduplication priority** — when the same underlying sequence shows up from both an MGnify and a third-party source, the MGnify record always wins, regardless of the `genome`/`metagenome` value on either side; the `genome` vs. `metagenome` preference only breaks ties _within_ a single origin (i.e. it decides between two MGnify records, or between two third-party records, but an MGnify `genome` record still outranks a third-party `metagenome` one). Biome and other metadata from every contributing record are merged into the final entry regardless of which one wins.

## Running without MGnify input

Third-party input does not depend on MGnify input at all — a run can be launched with `--third_party_input` alone and no `--input`. In that case, renaming starts its accession numbering from scratch (there is no MGnify range to continue from), and every downstream step (separation, quality evaluation, deduplication, clustering, annotation) proceeds exactly as described above, just without any MGnify-sourced records to merge against. Providing only `--input`, only `--third_party_input`, or both together are all valid, supported ways to run the pipeline.

## Run pipeline

```bash
nextflow run main.nf \
    -resume \
    -profile <appropriate profile> \
    -c <appropriate.config> \
    --outdir <OUTDIRNAME> \
    --third_party_input third_party_samplesheet.csv
```

To combine third-party data with an MGnify catalogue in the same run, add `--input` alongside `--third_party_input` (see [MGnify usage](mgnify_usage.md#run-pipeline) for the MGnify-specific options, such as `--rename_accession`, `--start_accession` and `--end_accession`):

```bash
nextflow run main.nf \
    -resume \
    -profile <appropriate profile> \
    -c <appropriate.config> \
    --outdir <OUTDIRNAME> \
    --input samplesheet.csv \
    --third_party_input third_party_samplesheet.csv
```

## Outputs

Pipeline results are written to the specified `OUTDIRNAME`, following the structure described in the [output documentation](output.md). Third-party sequences that fail FASTA validation are listed separately in `invalid_fastas.tsv` at the top of the output directory, rather than appearing anywhere in the catalogue outputs.
