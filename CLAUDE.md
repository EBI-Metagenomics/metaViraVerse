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
- `conf/modules.config` — all `ext.args`/`ext.prefix` and every `publishDir`; output paths and file names are defined here, not in modules.
- `development/` — local scratch data and HTML report prototypes; not part of the pipeline.

## Data flow and identifiers (needs several files to see)

1. **Rename** (`bin/rename_contigs.py`): every input sequence gets a catalogue ID `<rename_accession><n>` (`MGYV` + 10-digit zero-padded for MGnify). It writes the rename map with columns `original, temporary, short, biome, type, definition`: `temporary` = new ID, `short` = original ID up to the first space, `type` = `genome`/`metagenome`, `definition` = `virus`/`prophage`/`plasmid` (prefixed `third_party_` for third-party rows). Downstream code maps IDs back through this file.
2. **Separate** (`bin/separate_sequences.py`) into virus/prophage/plasmid using the map's `definition`. CheckV and barrnap run on viruses and prophages only.
3. **Choose/deduplicate** (`bin/choose_sequences.py`): deduplicates across all three categories by sequence SHA-256, with priority MGnify-metagenome > MGnify-genome > third-party-metagenome > third-party-genome; biomes of duplicates are merged. Low-quality records go to the excluded table.
4. **Viruses** = viral sequences + prophages pooled, clustered with vclust (95% ANI / 85% coverage; `conf/modules.config` `VCLUST_*`), then **only cluster representatives** are annotated: taxonomy (geNomad, VITAP, ViPhOGs), host (iPHoP; SpacePHARER with custom spacers), lifestyle (BACPHLIP), proteins (HMMER, AMR, PhaMMseqs). `BUILD_FINAL_GFF` merges annotations into one GFF. Protein steps live in `subworkflows/local/proteins_subwf` (`PROTEINS_PROCESSING`), which is called from both the viral and the plasmid subworkflow with different skip flags.
5. **Plasmids** clustered separately (vclust: pairs with ≥70% ANI / ≥50% coverage, Leiden at gANI 0.35), representatives annotated with AMR, MMseqs2 protein clustering, MOB-suite `mob_typer` and plaSquid (R scripts); `bin/annotate_plasmid_gff.py` merges these into the GFF; `bin/update_mobility_stats.py` combines plaSquid and MOB-suite mobility counts.

Original MGnify sequence names look like `MGYG000519004_12|viral_sequence-1:32972` (genome, contig, type, coordinates); the source genome is the part before the first `_`, which is how scripts join to `--catalogues_metadata` (lineage, sample, study).

## Conventions

- **Optional inputs/steps** (`--skip_*`, `--predict_host_from_custom_spacers`, `--annotation_db`): skipped steps emit `channel.empty()`. Before passing such a channel into a process, map to the file and add `.ifEmpty([])`; in the module build the flag conditionally (`def x_arg = x ? "--x ${x}" : ""`); in the script make the argument `required=False` and handle its absence. Otherwise the downstream process never runs.
- After `.join(..., remainder: true)` the `.map { ... }` closure must take exactly as many values as the joined tuple has; update it whenever a join is added.
- `FIND_CONCATENATE` (used to merge chunked outputs) names its output exactly `ext.prefix`, with no extension unless the prefix has one — `publishDir` patterns and `saveAs` must account for that.
- When a module's output file name changes, update the matching `pattern`/`saveAs` in `conf/modules.config`; a mismatch silently stops publishing.
- Intermediate files are published only with `--save_intermediates` (`enabled: params.save_intermediates` in `conf/modules.config`); keep `docs/output.md` in sync with which outputs are always published.
- Missing values in final metadata tables are written as the string `missing` (upstream steps may write `NA` / `not-provided`).
- Scripts must read both plain and gzip-compressed inputs; many intermediate files are `.gz`.
- Chunked steps (CheckV, iPHoP, MOB-suite, VITAP) concatenate chunk outputs, which can repeat header lines — parsers should skip repeated headers.
