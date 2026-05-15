# metaViraVerse — Protein Structure Subworkflow

## Overview

The `PROTEIN_STRUCTURE` subworkflow adds structural annotations to viral cluster
representative proteins identified by metaViraVerse. It is fully optional —
activated only when `.faa` sequences are provided in the samplesheet.

Existing pipeline users who do not provide `.faa` are completely unaffected.

---

## Workflow steps

| # | Step | Tool | Purpose |
|---|------|------|---------|
| 1 | Filter & length-check | custom Python | Remove seqs >1000 aa or >5% ambiguous residues |
| 2 | Structure prediction | ESMFold (Meta) | Predict 3D structures without MSA; pLDDT confidence per residue |
| 3 | Quality filter | inline | Flag mean pLDDT < 0.7; lower-confidence predictions retained but marked |
| 4 | Structural search — BFVD | Foldseek | Find viral structural relatives invisible to sequence BLAST |
| 5 | Structural search — PDB | Foldseek (optional) | Broader hits; host-mimicry and non-viral homologs |
| 6 | ECOD annotation | HHsearch + ECOD HHM DB | X/H/T/F-group domain classification from Zenodo HHM DB |
| 7 | SCOP cross-ref | inline with ECOD | SCOP mappings parsed from ECOD output |
| 8 | Annotation merge | custom Python | Extends existing stats TSV with 9 structural columns |
| 9 | Structural landscape | ProteinCartography | Interactive UMAP of viral structural space |

---

## Enabling the subworkflow

Add `.faa` to your samplesheet and pass `--run_protein_structure` at runtime:

```bash
nextflow run EBI-Metagenomics/metaViraVerse \
    --input samplesheet.csv \
    --run_protein_structure \
    --esmfold_mode api \
    --plddt_cutoff 0.7 \
    --outdir results/
```

### Samplesheet format (with `.faa`)

```csv
sample,fastq_1,fastq_2,faa
sample1,sample1_R1.fastq.gz,sample1_R2.fastq.gz,sample1_reps.faa
sample2,sample2_R1.fastq.gz,sample2_R2.fastq.gz,
```

The `.faa` column is optional per row. Samples without `.faa` proceed through
the standard pipeline unchanged; a warning is emitted for missing `.faa`.

---

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `run_protein_structure` | `false` | Enable this subworkflow |
| `esmfold_mode` | `'api'` | `'api'` (development) or `'local'` (production; requires ~16 GB VRAM) |
| `bfvd_db` | `'BFVD'` | Path to Foldseek BFVD DB, or `'BFVD'` to stream via bfvd.foldseek.com |
| `plddt_cutoff` | `0.7` | Minimum mean pLDDT for ECOD annotation and ProteinCartography |
| `foldseek_search_pdb` | `false` | Also search against PDB100 |

---

## Output files

```
results/
└── cluster_reps/
    ├── structures/
    │   └── <rep_id>.pdb
    ├── structure_confidence.tsv
    ├── bfvd_hits.tsv
    ├── ecod_annotations.tsv
    ├── scop_annotations.tsv
    ├── viral_sequences_reps_struct_stats.tsv
    └── protein_map/
        └── protein_map.html
```

### New columns in `viral_sequences_reps_struct_stats.tsv`

| Column | Description |
|--------|-------------|
| `struct_predicted` | yes / no |
| `mean_plddt` | 0–1; >0.7 = reliable |
| `struct_status` | success / failed / not_run |
| `bfvd_top_hit` | Best BFVD structural homolog |
| `bfvd_tm_score` | TM-score (>0.5 = same fold) |
| `bfvd_evalue` | E-value of best BFVD hit |
| `bfvd_pident` | Sequence identity of best BFVD hit |
| `ecod_xgroup` | ECOD X-group (possible distant homology) |
| `ecod_hgroup` | ECOD H-group (probable homology) |
| `ecod_fgroup` | ECOD family |

---

## Storage requirements

| Resource | Size | Strategy |
|----------|------|----------|
| ESMFold weights | ~2.5 GB | HuggingFace cache; downloads once at runtime; set `HF_HOME` to scratch |
| BFVD DB | 9.1 GB | Streamed via `bfvd.foldseek.com` by default; no local download required |
| ECOD HHM DB | ~3 GB | Fetched from Zenodo record 13993145 at runtime; cached in Nextflow `workDir` |
| ProteinCartography | ~500 MB | Snakemake conda environments; no heavy database |

---

## Caveats

1. **ESMFold public API** caps at ~400 aa and is rate-limited. Use `--esmfold_mode local`
   for batches >100 proteins or sequences >400 aa.
2. **pLDDT interpretation**: mean pLDDT < 0.7 indicates unreliable prediction; treat
   structural conclusions from such proteins as preliminary.
3. **BFVD version**: results depend on DB version. Default streams current BFVD;
   for reproducibility, pin a local DB copy and record the download date in methods.
4. **ProteinCartography**: all-vs-all TM-score comparison is O(N²). For >1000
   cluster reps, runtime is hours; ensure ample memory.
5. **ECOD HHM DB**: fetched from Zenodo record 13993145 on first run; cached thereafter.

---

## Citations

If you use this subworkflow, please cite:

- **ESMFold**: Lin et al. *Science* 2023 (doi:10.1126/science.ade2574)
- **Foldseek**: van Kempen et al. *Nature Methods* 2024 (doi:10.1038/s41592-023-02119-x)
- **BFVD**: Kim et al. *Nucleic Acids Research* 2024 (doi:10.1093/nar/gkae1119)
- **ECOD**: Cheng et al. *Nucleic Acids Research* 2014 (doi:10.1093/nar/gkt1248)
- **ProteinCartography**: Bigge et al. *Arcadia Science* 2024
