# Methods

This document describes, step by step, what ViraVerse does to a set of input sequences: how records from different sources are merged, how their quality is assessed, how redundancy is removed, and how the resulting viral and plasmid sequences are clustered and characterised. It is meant as a companion to the development notes kept alongside the pipeline code, expanded here into full prose so that each stage of processing is explained together with the tool and rationale behind it.

## Input

The pipeline accepts two independent kinds of input, which are merged during pre-processing: **MGnify pre-processed data** and **third-party data**. The former is the output of the standard MGnify assembly/MAG analysis pipelines (viral and plasmid predictions already made with tools such as geNomad or VIRify), while the latter allows sequences from any external source (e.g. a public database dump, a collaborator's dataset, or a previously published catalogue) to be folded into the same run. See [MGnify docs](mgnify_usage.md) and [Third party docs](third_party_usage.md) for details on how to prepare a samplesheet for each.

MGnify input format:

| Column   | Required | Description                                                                                                                                                                                                                                                                         |
| -------- | -------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `id`     | Yes      | Unique identifier for the sample/assembly. We recommend using the ERZ accession when the MAG or assembly originates from ENA, so that results can always be traced back to the original submission.                                                                                 |
| `fna`    | Yes      | FASTA file with the nucleotide sequences for the regions listed in the GFF (i.e. the candidate viral/plasmid contigs, not the whole assembly).                                                                                                                                      |
| `gff`    | Yes      | GFF file describing the viral and plasmid records predicted for this sample. It may also contain CDS features for the selected regions, which are carried through the pipeline alongside the nucleotide sequences.                                                                  |
| `faa`    | No       | FASTA file with the protein sequences translated from the CDS regions in the GFF. When supplied, it is used downstream instead of re-predicting genes.                                                                                                                              |
| `source` | Yes      | `genome` (the sequence comes from a MAG) or `metagenome` (the sequence comes from an assembly). This label is later used to decide which copy of a duplicated sequence to keep — the same column, and the same meaning, as the third-party samplesheet's own `source` column below. |
| `biome`  | No       | Free-text metadata describing the sequence's environmental origin (for example: marine, soil, human gut). Retained purely as metadata and merged across duplicates.                                                                                                                 |

Third party input format:

| Column             | Required | Description                                                                                                                                                                                                                                                |
| ------------------ | -------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `id`               | Yes      | Unique identifier for the record.                                                                                                                                                                                                                          |
| `fna`              | Yes      | FASTA file with the nucleotide sequence(s) for the chosen `type`.                                                                                                                                                                                          |
| `type`             | Yes      | `virus`, `prophage`, or `plasmid`.                                                                                                                                                                                                                         |
| `source`           | No       | `genome` (the sequence comes from a MAG or isolate) or `metagenome` (the sequence comes from an assembly) — the same distinction MGnify's own `type` column makes, and used the same way downstream to decide which copy of a duplicated sequence to keep. |
| `biome`            | No       | Metadata describing the sequence's environmental origin (for example: marine, soil).                                                                                                                                                                       |
| `study_accession`  | No       | INSDC study accession (primary or secondary) associated with the data, kept as metadata.                                                                                                                                                                   |
| `sample_accession` | No       | INSDC sample accession associated with the data, kept as metadata.                                                                                                                                                                                         |

Third-party records go through their own preparation subworkflow (annotation, gene calling where needed) before being merged with the MGnify samplesheet, so that from the pre-processing stage onward both sources are treated uniformly. A run may supply only MGnify input, only third-party input, or both together.

## Pre-processing

Pre-processing is the stage where all inputs — regardless of whether they came from MGnify or from a third-party samplesheet — are brought together, split into viruses, prophages and plasmids, quality-checked, and de-duplicated (across all three categories) into one coherent collection of sequences ready for downstream analysis.

### Rename sequences

Before anything else, every incoming record is given a short, unique, pipeline-internal identifier. This matters because contig/MAG headers coming from different assemblies or external sources are not guaranteed to be unique (or even safely parsable), and several downstream tools are sensitive to long or unusual FASTA headers. This step renames every sequence to `{prefix}N` (for example `seq1`, `seq2`, …), also updating the matching `ID` attributes inside the GFF so that features stay linked to the right contig.

- MGnify sequences are renamed first; third-party sequences are renamed immediately afterwards, continuing the very same accession counter rather than starting over — so the whole catalogue ends up with one contiguous numbering range regardless of how many sources contributed to it. This works even when a run supplies only one of the two sources (MGnify-only or third-party-only).
- The prefix defaults to `seq` and is controlled by `--rename_accession`; for MGnify catalogue generation it is recommended to set this to `MGYV` so identifiers match the MGnify accession convention.
- The numeric part of the accession can be bounded with `--start_accession` and `--end_accession`, which is useful when a catalogue is built incrementally and new sequences need to continue numbering from where a previous batch left off.
- A mapping file (`combined.tsv`) is produced alongside the renamed, combined FASTA/GFF, recording, for every new temporary name: the original sequence name, its `biome`, its `type`, and its `definition`. This mapping is consulted by nearly every later step (separation, quality filtering, deduplication, cluster extraction) whenever the original identity or provenance of a sequence needs to be recovered.
  - `type` always means `genome` (MAG/isolate) or `metagenome` (assembly), regardless of where the sequence came from: both the MGnify and the third-party samplesheets carry this value in their own `source` column, and it is copied verbatim into the mapping file's `type` column for either origin. Keeping this single, consistent vocabulary is what lets the deduplication step (below) compare an MGnify and a third-party record on equal footing.
  - `definition` records the sequence's biological category — `virus`, `prophage`, or `plasmid`. For MGnify records this is derived automatically by testing the sequence's original name against the `--viral_sequence_identifier`, `--prophage_identifier` and `--plasmid_identifier` patterns (the same three patterns used to physically separate sequences afterwards). For third-party records it is simply copied from the samplesheet's own `type` column, but prefixed with `third_party_` (e.g. `third_party_prophage`) — that prefix is what lets later steps recognise a record as third-party-derived even though its `type` (genome/metagenome) column looks identical in shape to an MGnify record's.

### Separate into viruses, prophages and plasmids

Immediately after renaming, the combined FASTA, GFF and protein (FAA) files are split into three groups — viruses, prophages, and plasmids — using the `definition` column recorded in the mapping file during renaming. A sequence whose definition is `third_party_prophage` is filed into the same "prophage" group as one whose definition is plainly `prophage`; the `third_party_` prefix is stripped for the purpose of bucketing, so that MGnify and third-party sequences of the same biological category are always treated together from this point on.

- All three file types are split together: the GFF's CDS features are traced back to their parent sequence so that, for every sequence kept in a given group, its own GFF records and its own predicted proteins travel with it into that group's `viruses.fna`/`.gff`/`.faa` (and likewise for `prophages.*` and `plasmids.*`).
- Splitting this early — before quality evaluation and deduplication, rather than after — means every downstream step can be scoped to just the categories it actually applies to (for example, plasmids skip CheckV and Barrnap entirely, see below), and the mapping-derived category travels alongside the sequence data itself rather than needing to be re-derived later.

### Quality evaluation

Sequence quality is assessed for viral and prophage sequences (plasmids are not run through CheckV) using [CheckV](https://bitbucket.org/berkeleylab/checkv/src/master/), which estimates completeness and contamination for viral and proviral genomes by comparing them against a database of reference viral genomes and a database of closely related proteins.

- The virus and prophage FASTA files are first combined and then split into chunks (size controlled by `--nucleotide_fasta_chunksize`, 25,000 sequences by default) so that CheckV's `end_to_end` module can be parallelised across many small jobs rather than run once on a potentially very large file.
- Each chunk is run through `checkv end_to_end`, which internally executes CheckV's three modules — contamination removal, completeness estimation (AAI- or HMM-based), and complete-genome identification — against `--checkv_db`.
- The per-chunk `quality_summary.tsv` outputs are concatenated back into a single table describing, per sequence: whether it is a provirus and its proviral length, gene and viral-gene counts, CheckV/MIUVIG quality tier, completeness (and the method used to estimate it), contamination percentage, k-mer frequency, and any warnings raised by CheckV.
- This table becomes one of the primary inputs to the deduplication/filtering step below, where sequences flagged as low-confidence are removed from the catalogue.

### rRNA prediction

To avoid keeping ribosomal/host-derived contamination in the viral catalogue, the same combined virus+prophage sequence set (plasmids are excluded here too) is screened for ribosomal and transfer RNA genes using [Barrnap](https://github.com/tseemann/barrnap), which uses HMMER profiles built from bacterial, archaeal and eukaryotic rRNA models to locate rRNA genes in nucleotide sequences.

- Barrnap is run against the bacterial rRNA model, since the sequences being screened are expected to derive from bacterial (or archaeal-adjacent) hosts rather than being rRNA genes themselves.
- Any `rRNA`, `tRNA` or `tmRNA` feature reported in the resulting GFF flags its parent sequence as rRNA/tRNA-bearing.
- This screening step can be skipped entirely with `--skip_rrna_detection` (for example on a very well curated input set, or to save runtime), in which case downstream steps treat rRNA status as "not provided" rather than filtering on it.
- A sequence flagged as carrying rRNA/tRNA/tmRNA is excluded from the catalogue in the deduplication/filtering step below — plasmids were never screened in the first place, so this exclusion only ever applies to virus/prophage sequences.

### Sequences deduplication

The same underlying sequence is often predicted independently in more than one sample, recovered both as part of an assembly and as part of a MAG derived from that assembly, or even classified differently by different predictions (for example called a free virus in one sample and a prophage in another). This step compares sequences **across all three categories at once**, collapses duplicates into one representative record per unique sequence, removes low-confidence calls, and produces the metadata table and filtered FASTA/GFF/FAA that everything downstream is built on.

- Every sequence, from every one of the three category inputs, is hashed (SHA256, on the uppercased nucleotide string). Two records that hash identically are treated as the same underlying sequence regardless of which sample, source, or even category (virus/prophage/plasmid) they were originally filed under — this is what catches a sequence that was called a prophage in one prediction and a free virus in another.
- When a sequence appears more than once, only one copy is retained, chosen by a two-level priority: an MGnify-derived record always outranks a third-party one, and within either origin a `metagenome` (assembly) record outranks a `genome` (MAG) one. In order, that gives: MGnify+metagenome, MGnify+genome, third-party+metagenome, third-party+genome. Biome labels from every contributing source are merged and kept regardless of which record wins, and the deduplicated sequence's final category (virus/prophage/plasmid) becomes whichever category the _winning_ record came from — so a sequence can end up re-classified from how it first appeared if a higher-priority duplicate placed it in a different category.
- Two extra pieces of information are attached to every retained sequence from the earlier steps: whether it carries an rRNA/tRNA/tmRNA annotation (from Barrnap), and its CheckV quality metrics — both apply only to virus/prophage sequences, since plasmids were never screened.
- Filtering then removes sequences that are unlikely to be genuine viral calls:
  - Any non-plasmid sequence carrying an rRNA/tRNA/tmRNA annotation is excluded.
  - Any non-plasmid sequence whose CheckV quality is `Not-determined`, with zero viral genes detected and a k-mer frequency ≥ 1.0 (a combination that is characteristic of host contamination rather than a genuine, just-hard-to-classify virus), is excluded.
- Because GFF and protein records are carried alongside their sequence throughout, dropping a sequence — whether because it lost a duplicate comparison or failed the quality filter — also drops its corresponding GFF features and predicted proteins from the output, rather than leaving orphaned annotations behind.
- The step produces: a full metadata table for every unique sequence seen (`*_metadata.tsv`, now including a `definition` column recording the final virus/prophage/plasmid/third-party status), the filtered FASTA/GFF/FAA/TSV kept for downstream processing — both per-category (`*_virus_filtered.*`, `*_prophage_filtered.*`, `*_plasmid_filtered.*`) and combined across all three (`*_filtered.*`) — and a table of excluded sequences together with the reason each was dropped (`*_excluded.tsv`), kept so that filtering decisions remain auditable. Free viral sequences and prophages are merged again downstream and processed together as "viruses", while plasmids follow a separate branch.

## Sequences processing: Plasmids

Once plasmid sequences have been separated out, they are clustered into representative genomes and then characterised with tools specialised for plasmid biology (replicon/mobility typing and plasmid-specific gene detection).

### Clustering plasmids

Plasmid sequences are clustered into representative genomes ("plasmid operational taxonomic units", POTUs) using [vclust](https://github.com/refresh-bio/vclust/wiki/6-Use-cases#68-cluster-plasmid-genomes-into-potus), following the use case described in the vclust wiki for grouping plasmid genomes.

- Clustering is based on global ANI (gANI), computed by `vclust align` after `vclust prefilter` has discarded genome pairs too dissimilar to be worth aligning, and then thresholded by `vclust cluster`.
- The similarity threshold used for plasmids is a minimum gANI of 0.35 — considerably looser than the threshold used for viruses (see below), reflecting the greater sequence diversity typically observed within a single plasmid lineage.
- If `--cluster_vclust` is disabled, an alternative, simpler pipeline is used instead: an all-vs-all BLASTN search (via `makeblastdb`/`blastn`), followed by pairwise ANI calculation (`anicalc`) and a UCLUST-style greedy centroid clustering (`aniclust`), both adapted from CheckV's companion scripts.
- Whichever method is used, the output is a table assigning every plasmid sequence to a cluster and identifying one representative per cluster; only these representatives are carried forward into the plasmid-specific analyses below (as well as into the website indexing step).

### MOB

Representative plasmid sequences are typed with [MOB-suite's `mob_typer`](https://github.com/phac-nml/mob-suite#mob-typer), which provides in silico predictions of the replicon family, relaxase type, mate-pair formation (MPF) type, and predicted transferability of each plasmid.

- Using a combination of sequence biomarkers and MOB-cluster codes, `mob_typer` also estimates the observed host range of a plasmid from its replicon, relaxase, and cluster assignment, combined with information mined from the literature, to predict the taxonomic rank at which the plasmid is likely to be stably maintained. It does **not**, however, attempt source (i.e. specific host organism) attribution.
- This gives the catalogue, for every representative plasmid, a mobility classification (conjugative, mobilisable, or non-mobilisable) alongside its predicted replicon and relaxase types, which is valuable metadata for interpreting how the plasmid may spread between hosts.

### plaSquid

[plaSquid](https://github.com/mgimenez720/plaSquid) is run to complement MOB-suite with additional plasmid-specific gene detection. The original plaSquid pipeline was re-written and re-integrated into ViraVerse as individual Nextflow modules, covering the replicon search, plasmid-associated RNA search, incompatibility-group search and relaxase/mobilisation search modes (`repsearch`, `rnasearch`, `incsearch` and `mobsearch`).

- `repsearch` screens the predicted proteins against a replicon domain database (`--repsearch_db`, filtered with `--repfilter_db`) to detect replication-initiation genes characteristic of plasmids.
- `rnasearch` screens the nucleotide sequence itself (`--rna_inc_db`) for plasmid-associated RNA features (for example RNAs involved in copy-number control).
- `incsearch` combines the protein search results with the RNA search candidates to classify plasmids into incompatibility groups (`--incsearch_db`).
- `mobsearch` screens the predicted proteins for mobilisation/relaxase genes (`--mobsearch_db`).
- The results of all four searches are finally combined into a single per-plasmid classification table, giving a second, complementary line of evidence to the MOB-suite predictions above.

## Sequences processing: Viruses

Viral sequences and prophages, once merged back together, are clustered into representative viral genomes and then run through host, lifestyle and taxonomy prediction before the final annotated catalogue is assembled.

### Clustering viruses and prophages

Viral (including prophage) sequences are dereplicated into representative genomes using [vclust](https://github.com/refresh-bio/vclust/wiki/6-Use-cases#63-dereplicate-viral-contigs-into-representative-genomes), following the standard use case for building viral operational taxonomic units (vOTUs).

- Clustering here is based on ANI rather than gANI, with a minimum identity of 0.95 (95%) between a sequence and its cluster representative — the community-standard threshold (in line with MIUVIG recommendations) for grouping viral genomes into species-level vOTUs.
- As with plasmid clustering, the same `--cluster_vclust` switch controls whether vclust or the BLASTN + `anicalc`/`aniclust` fallback pipeline is used.
- Only cluster representatives are propagated to the rest of the viral analysis branch (host detection, lifestyle prediction, taxonomy assignment, protein annotation), keeping the compute-intensive downstream steps focused on non-redundant genomes while every member sequence remains recorded in the clustering table.

### Lifestyle

The lifestyle of each representative phage/virus genome — whether it is predicted to be **virulent** (strictly lytic) or **temperate** (capable of lysogeny) — is predicted with [BACPHLIP](https://github.com/adamhockenberry/bacphlip), which classifies bacteriophage genomes from the presence/absence pattern of a curated set of conserved protein domains associated with the lysogenic cycle (e.g. integrases, excisionases, repressors).

- Representative sequences are chunked (using the same `--nucleotide_fasta_chunksize` setting as quality evaluation) before being run through BACPHLIP, so the prediction scales to large catalogues.
- Per-chunk results are concatenated into a single lifestyle table, which is later merged into the final annotated GFF alongside protein annotations, so that lifestyle can be inspected side-by-side with gene content for every representative genome.
- BACPHLIP runs inside a dedicated container, built from the [Dockerfile](../container/bacphlip/Dockerfile) included in this repository, pinned to BACPHLIP v0.9.3-alpha.

### Taxonomy assignment (DNA sequence based)

Two independent, optional tools assign taxonomy to representative viral/prophage genomes directly from their nucleotide sequence, without needing predicted proteins first:

- **[VITAP](https://github.com/DrKaiZheng/VITAP)** (skippable with `--skip_vitap`) classifies the chunked representative nucleotide sequences against `--vitap_db`, a high-precision meta-omic viral classification model, and its single best lineage call per sequence is kept.
- **[geNomad](https://github.com/apcamargo/genomad)** (skippable with `--skip_genomad`) classifies each representative sequence from marker-gene composition and a deep neural-network model trained on curated viral genomes, assigning both a virus/plasmid/chromosome call and, for viral sequences, a taxonomic lineage.
- VITAP and geNomad are run independently and are **not** reconciled into a single consensus call: each tool's raw output is instead combined with the deduplicated sequence metadata into its own taxonomy table, so that the two lineage assignments can be compared side by side (together with the protein-based ViPhOGs assignment described below) rather than one silently overriding the other.

### Taxonomy visualisation

For every taxonomy source (VITAP, geNomad and ViPhOGs), a set of complementary visual summaries is generated from the combined taxonomy-and-metadata table:

- **[Krona](https://github.com/marbl/Krona)** (`KRONA_KTIMPORTTEXT`) renders the taxonomic composition of the catalogue as an interactive, zoomable radial chart, letting a reader drill down from domain/realm level to species-level vOTUs.
- A **Sankey plot** renders the same lineage assignments as a flow diagram, which is often more readable than a Krona chart for tracing how sequences split across a small number of intermediate taxonomic ranks.
- **[iTOL](https://itol.embl.de/)**-ready tree annotation files are generated so that a phylogenetic or clustering tree of the catalogue can be uploaded to iTOL and coloured/annotated by taxonomy directly in the browser.

### Host detection

Two complementary, optional approaches link representative viral/prophage genomes to their likely prokaryotic hosts:

- **[iPHoP](https://bitbucket.org/srouxjgi/iphop)** (skippable with `--skip_iphop`) integrates several existing host-prediction signals (including CRISPR spacer matches, k-mer composition similarity, and prophage sequence similarity to reference genomes) behind a single machine-learning framework, and reports a genus-level and a genome-level host call with a confidence score for each representative sequence. Representatives are chunked with their own, smaller chunk size (`--nucleotide_fasta_chunksize_iphop`, 5,000 sequences by default, reflecting iPHoP's higher per-sequence cost) and run against `--iphop_db`; per-chunk genus- and genome-level tables are concatenated back together.
- **[SpacePHARER](https://github.com/soedinglab/spacepharer)**, run only when `--predict_host_from_custom_spacers` is enabled, searches representative genomes directly against a user-supplied, pre-generated database of host CRISPR spacers (`--custom_spacers_fasta` / `--custom_spacers_metadata`). This is useful when a project has its own curated spacer set (for example spacers extracted from paired isolate genomes) that is more relevant to the sampled environment than iPHoP's general-purpose training data. SpacePHARER compares spacers and phage sequences at the protein level for extra sensitivity to these short, fast-evolving sequences, and matches are combined with the supplied spacer metadata into a single host-prediction table.

## Proteins processing: Viruses

### Protein calling for Third party data

Unlike MGnify input (which always arrives with predicted proteins), third-party records may only supply a nucleotide sequence. Before such records can go through protein-based annotation (ViPhOGs taxonomy, AMR detection, functional annotation), genes are predicted with **[Pyrodigal](https://github.com/althonos/pyrodigal)**, a Python/Cython reimplementation of Prodigal. Input FASTA files are first decompressed and validated with [FALint](https://github.com/GallVp/fa-lint) (files that fail validation are logged to `invalid_fastas.tsv` and excluded from the run, rather than failing it outright), then the survivors are run through Pyrodigal to produce a GFF of predicted CDS features and the corresponding translated protein FASTA, so that from this point onward third-party and MGnify records carry the same `gff`/`faa` pair. See [Third-party data usage](third_party_usage.md) for the full samplesheet format and how these records are then merged with the MGnify catalogue.

### Taxonomy assignment (Aminoacid sequence based)

**ViPhOGs** (skippable with `--skip_viphogs`), a set of profile HMMs originally developed for and reused from the [VIRify](https://github.com/EBI-Metagenomics/emg-viral-pipeline) pipeline, is searched against representative predicted proteins (`--viphog_db`, with `--additional_model_data` and `--ncbi_db` supplying supporting reference data, and an optional `--factor_file` for taxon-weighting). A per-genome taxonomic assignment is derived from the pattern of ViPhOG hits along each representative genome, giving a protein-based lineage call that complements the two nucleotide-based calls (VITAP, geNomad) described above.

### AMR annotation

Antimicrobial-resistance genes are detected with up to three independent, optional tools, run over the same representative protein set: [AMRFinderPlus](https://github.com/ncbi/amr) (curated reference-gene search against the NCBI Reference Gene Catalog; disabled by default — enable with `--skip_amrfinderplus false`), [DeepARG](https://github.com/gaarangoa/deeparg) (deep-learning classification of resistance genes, `--skip_deeparg`) and [RGI](https://github.com/arpcard/rgi) against the CARD database (`--skip_rgi`). Each tool's raw output is standardised with hAMRonization and the (up to three) sets of calls are integrated into a single AMR GFF track, so that a gene flagged by more than one detector is easy to spot. See the [AMR annotation subworkflow docs](https://ebi-metagenomics.github.io/nf-modules/subworkflows/ebi-metagenomics/amr_annotation/) for the full module-by-module breakdown.

### Hmmsearch annotation

General protein function is annotated by searching every representative protein against a set of curated profile-HMM databases with **[HMMER](http://hmmer.org/)** `hmmsearch`. The idea, and the bundled HMM databases themselves, are taken from the [MetaCerberus](https://github.com/raw-lab/metacerberus) tool — ViraVerse simply runs `hmmsearch` directly against these (in some cases updated) databases rather than depending on the full MetaCerberus pipeline. Proteins are chunked beforehand (`--protein_annotation_fasta_chunksize`, 50,000 sequences by default) so that the search parallelises across many HMM databases and protein chunks at once; per-database hit tables are concatenated and combined with each database's descriptive metadata into a single functional-annotation summary, which is later merged into the final per-genome GFF alongside the lifestyle and AMR tracks.

### Protein clustering

As an optional step (`--phammseqs`, off by default), representative proteins are grouped into "phamilies" of related sequences with **[PHAMMSeqs](https://github.com/chg60/phammseqs)**, an MMseqs2-based clustering wrapper originally developed for bacteriophage comparative genomics. This is mainly useful for downstream comparative-genomics work rather than for the core catalogue itself, so it is left disabled unless specifically requested.

## Protein structures

Structural prediction/annotation for representative proteins is still in development and not yet part of a released version of the pipeline; see [PR #19](https://github.com/EBI-Metagenomics/metaViraVerse/pull/19) for the current state of that work.
