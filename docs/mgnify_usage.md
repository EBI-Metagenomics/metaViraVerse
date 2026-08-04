# MGnify usage

## Data preparation

Viruses and plasmids for the viral catalogue are detected using [emg-viral-pipeline](https://github.com/EBI-Metagenomics/emg-viral-pipeline) (VIRify) and [mobilome-annotation-pipeline](https://github.com/EBI-Metagenomics/mobilome-annotation-pipeline) (MAP), run as part of:

- [MAG catalogue](https://www.ebi.ac.uk/metagenomics/browse/genomes) generation
- [Assembly analysis pipeline](https://github.com/EBI-Metagenomics/assembly-analysis-pipeline) (ASA) runs

**Important:**

VIRify+MAP use pre-defined names for viruses (**viral_sequence**), prophages (**prophage**) and plasmids (**plasmid**). Those names are used in 2 places:

- Fetch [script](../scripts/collect_data_from_catalogues.py) to grep sequences of interest from MAG's GFF file.
- pipeline [separate_sequences](../bin/separate_sequences.py) step to split sequences into groups.

Make sure those names haven't changed (check MAP). Pipeline execution can regulate those names via `--viral_sequence_identifier`, `--prophage_identifier` and `--plasmid_identifier`.

> [!NOTE]
> Viral catalogue generation pipeline has only been tested on data from MAG catalogues so far.

### Fetch data from MGnify genomes

**Mandatory step:**

- Collect data from existing catalogues (this also fetches each catalogue's metadata and generates the pipeline input samplesheet)

**Optional step:**

- Collect CRISPR spacers

#### Collect data from existing catalogues

Make sure the catalogues you use included both VIRify and MAP in their generation. List of catalogues available as of 21/07/2026:

- Barley rhizosphere
- Maize rhizosphere
- Tomato rhizosphere
- Soil
- Marine
- Marine sediment
- Human skin
- Human vaginal
- Chicken gut
- Non-model fish gut
- Mouse gut
- Honey bee gut
- Sheep rumen

Choose the catalogues you want to use and find their locations on `/nfs/public/`. The last two components of each path (catalogue name and version, e.g. `human-vaginal/v1.0`) must match the corresponding catalogue location on the [MGnify genomes FTP](https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/), since that's where metadata is fetched from.

Then run the fetching script [`collect_data_from_catalogues.py`](../scripts/collect_data_from_catalogues.py) (make sure you're in the correct queue to access NFS):

```commandline
usage: collect_data_from_catalogues.py [-h] -p CATALOGUE_PATH [CATALOGUE_PATH ...] -o OUTPUT_PATH

Script searches for viral records in catalogue(s) GFFs and greps corresponding nucleotide and protein sequences.

options:
  -h, --help            show this help message and exit
  -p, --catalogue-path CATALOGUE_PATH [CATALOGUE_PATH ...]
                        Path(s) to NFS location of catalogue
  -o, --output-path OUTPUT_PATH
                        Path to save results (filtered gff, fna, faa)

        Script takes as input path(s) to catalogue(s) and, for each one, creates a filtered
        GFF/FNA/FAA and downloads the catalogue's metadata from the MGnify genomes FTP.
        Once all catalogues are processed, it also writes a combined metadata TSV and a
        pipeline input samplesheet covering every catalogue.
```

Example:

```bash
python3 collect_data_from_catalogues.py \
  -p nfs/catalogue_1/v1.0 nfs/catalogue_2/v1.0 \
  -o results
```

<details>
<summary>Outputs description</summary>

- `<CATALOGUE_NAME>_<VERSION>_viral.fna` — FASTA file with viral sequences (viral sequences, prophages, plasmids) found in `CATALOGUE_NAME_VERSION`.

  Header format example:

  ```
  >MGYG000517142_2 prophage|311953:328572
  ```

  Where:
  - `MGYG000517142` is the MGnify genome identifier
  - `MGYG000517142_2` is the contig identifier
  - `prophage` is the sequence type
  - `311953:328572` is the region coordinates

- `<CATALOGUE_NAME>_<VERSION>_viral.faa` — FASTA file with proteins detected within the chosen viral regions.

  Header format example:

  ```
  >MGYG000517142_00528 Alpha-xylosidase
  ```

  Where:
  - `MGYG000517142_00528` is the protein identifier
  - `Alpha-xylosidase` is the protein name

- `<CATALOGUE_NAME>_<VERSION>_viral.gff` — GFF with all viral records and the CDS records detected within those regions:

  ```
   MGYG000517142_2 geNomad prophage        311953  328572  .       .       .       ID=MGYG000517142_2|prophage-311953:328572;mobile_element_type=prophage;taxonomy=Viruses%3BDuplodnaviria%3BHeunggongvirae%3BUroviricota%3BCaudoviricetes%3B%3B
   MGYG000517142_2 Prodigal:002006 CDS     18      1355    .       +       0       ID=MGYG000517142_00528;eC_number=3.2.1.177;Name=yicI;Dbxref=COG:COG1501;gene=yicI;inference=ab initio prediction:Prodigal:002006,similar to AA sequence:UniProtKB:P31434;locus_tag=MGYG000517142_00528;product=Alpha-xylosidase;eggNOG=1194165.CAJF01000023_gene3185;cog=G;kegg=ko:K01811;pfam=PF01055,PF21365;interpro=IPR000322,IPR013780,IPR017853,IPR048395,IPR050985
   MGYG000517142_2 Prodigal:002006 CDS     1352    2776    .       +       0       ID=MGYG000517142_00529;eC_number=3.2.1.21;inference=ab initio prediction:Prodigal:002006,similar to AA sequence:UniProtKB:P94248;locus_tag=MGYG000517142_00529;product=Bifunctional beta-D-glucosidase/beta-D-fucosidase;eggNOG=1194165.CAJF01000023_gene3186;cog=G;kegg=ko:K05350;pfam=PF00232;interpro=IPR001360,IPR017736,IPR017853
  ```

  Each record follows the structure: viral region, followed by its predicted CDS entries.

  **Important:**
  - **Region IDs** correspond to whole sequence headers in `<CATALOGUE_NAME>_<VERSION>_viral.fna`.
  - **CDS IDs** correspond to protein identifiers in `<CATALOGUE_NAME>_<VERSION>_viral.faa`.

- `<CATALOGUE_NAME>_<VERSION>_metadata.tsv` — per-catalogue metadata, downloaded from `https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/<CATALOGUE_NAME>/<VERSION>/genomes-all_metadata.tsv`.

- `genomes-all_metadata.tsv` — all per-catalogue metadata files concatenated, keeping a single header row. This is the file to pass to `--catalogues_metadata` when running the pipeline.

- `samplesheet.csv` — pipeline input samplesheet with one row per catalogue, following [`assets/schema_input.json`](../assets/schema_input.json). See [Input samplesheet](#input-samplesheet) below.

</details>

#### Collect CRISPR spacers (optional)

One of the pipeline's steps is viral host detection. To find a possible host among MGnify genomes, you need to collect all spacers predicted by CRISPRCasFinder. This can be done with the [`collect_crispr_spacers_from_catalogues.py`](../scripts/collect_crispr_spacers_from_catalogues.py) script (make sure you're in the correct queue to access NFS):

```commandline
usage: collect_crispr_spacers_from_catalogues.py [-h] -p CATALOGUE_PATH [CATALOGUE_PATH ...] -o OUTPUT_PATH --prefix PREFIX

Script searches for CRISPRCasFinder results in catalogue(s) and greps CRISPR spacer information.

options:
  -h, --help            show this help message and exit
  -p, --catalogue-path CATALOGUE_PATH [CATALOGUE_PATH ...]
                        Path to NFS location of catalogue(s)
  -o, --output-path OUTPUT_PATH
                        Path to save results
  --prefix PREFIX       Output filename
```

Example:

```bash
python3 collect_crispr_spacers_from_catalogues.py \
  -p nfs/catalogue_1/v1.0 nfs/catalogue_2/v1.0 \
  -o results_spacers \
  --prefix all_catalogues
```

<details>
<summary>Outputs description</summary>

- `PREFIX_crispr.fasta` — multi-FASTA file containing spacer sequences.
- `PREFIX_crispr.tsv` — TSV with spacer metadata: ID, name, and parent (genome, protein ID, and fragment coordinates).

Example:

```commandline
crispr_id       name    parent
sp_19345        spacer_19345_29 MGYG000518600_13_19319_19399
```

</details>

## Input samplesheet

The "Collect data from existing catalogues" step generates `samplesheet.csv` in the output directory automatically — no manual step needed. It contains one row per catalogue, following [`assets/schema_input.json`](../assets/schema_input.json):

```csv
id,gff,fna,faa,type,biome
human-vaginal_v1.0,/path/to/results/human-vaginal_v1.0_viral.gff,/path/to/results/human-vaginal_v1.0_viral.fna,/path/to/results/human-vaginal_v1.0_viral.faa,genome,human-vaginal
```

| Column  | Required | Description                                                                                                  |
| ------- | -------- | ------------------------------------------------------------------------------------------------------------ |
| `id`    | Yes      | Unique identifier. We recommend using the ERZ accession if the MAG or assembly originates from ENA.          |
| `gff`   | Yes      | GFF file containing viral and plasmid records. May also contain CDS records for the selected regions.        |
| `fna`   | Yes      | FASTA file with nucleotide sequences for the selected regions in the GFF.                                    |
| `faa`   | No       | FASTA file with protein sequences for the CDS regions in the GFF.                                            |
| `type`  | Yes      | `genome` (for a MAG source) or `metagenome` (for an assembly source), describing the origin of the sequence. |
| `biome` | No       | Metadata describing the sequence's environmental origin (for example: marine, soil).                         |

## Run pipeline

```bash
nextflow run main.nf \
    -resume \
    -profile <appropriate profile> \
    -c <appropriate.config> \
    --outdir <OUTDIRNAME> \
    --input samplesheet.csv \
    --catalogues_metadata genomes-all_metadata.tsv \
    --rename_accession MGYV \
    --start_accession <first identifier number, e.g. 30, used when renaming into MGYV space, e.g. MGYV0000030> \
    --end_accession <last identifier number, e.g. 40, used when renaming into MGYV space, e.g. MGYV0000040> \
    --predict_host_from_custom_spacers <if CRISPR spacers were fetched in the "Collect CRISPR spacers" step> \
    --custom_spacers_fasta PREFIX_crispr.fasta <if CRISPR spacers were fetched in the "Collect CRISPR spacers" step> \
    --custom_spacers_metadata PREFIX_crispr.tsv <if CRISPR spacers were fetched in the "Collect CRISPR spacers" step> \
    --phammseqs <if you want to run protein clustering>
```

## Outputs

Pipeline results are written to the specified `OUTDIRNAME`, following the structure described in the [output documentation](output.md).
