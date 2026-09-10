## Databases for Viruses and Prophages

### Download mandatory databases

| Database               | Notes                                                                                                                                                                                                                                                                                                                     | Input argument         | Link                                                           | Version          | Publication                                                            |
|------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|------------------------|----------------------------------------------------------------|------------------|------------------------------------------------------------------------|
| CheckV                 | CheckV database                                                                                                                                                                                                                                                                                                           | `--checkv_db`          | https://portal.nersc.gov/CheckV/checkv-db-v1.5.tar.gz          | 1.5              | https://www.nature.com/articles/s41587-020-00774-7                     |
| Protein annotation DBs | Databases were originally taken from MetaCerberus tool and updated to latest available version. <br/>All databases stored in format `<DB>.hmm.gz` and `<DB>.tsv`. <br/>`TSV` file contains columns with metadata, for example, `ID` and `Function`. All files for all annotation databases should be stored in one folder | `--annotation_db`      | https://github.com/raw-lab/MetaCerberus#metacerberus-databases | 1.4.0            | https://academic.oup.com/bioinformatics/article/40/3/btae119/7616988   |
| VOG                    | The Virus Orthologous Groups Database (VOGDB) is a multi-layer database that progressively groups viral genes into groups connected by increasingly remote homology.                                                                                                                                                      | part of annotation dbs | https://vogdb.org/                                             | 236 (2026-02-25) | https://doi.org/10.3390/v16081191                                      |
| PHROG                  | Prokaryotic Virus Remote Homologous Groups database                                                                                                                                                                                                                                                                       | part of annotation dbs | https://phrogs.lmge.uca.fr/                                    | 2022-06-15       | https://doi.org/10.1093/nargab/lqab067                                 |
| PGAPfams               | PGAPfams are a collection of protein family Hidden Markov Models (HMMs) used by the NCBI Prokaryotic Genome Annotation Pipeline (PGAP) to automatically identify and assign functional annotations to microbial protein sequences                                                                                         | part of annotation dbs | https://ftp.ncbi.nlm.nih.gov/hmm/current/                      | 2026-01-27       | https://doi.org/10.1093/nar/gkw569                                     |
| NFixDB                 | Nitrogen Fixation DataBase                                                                                                                                                                                                                                                                                                | part of annotation dbs | https://github.com/raw-lab/NFixDB                              | 2024-01-22       | https://academic.oup.com/nargab/article/6/2/lqae063/7688810?login=true |
| GVDB                   | Giant Virus Database                                                                                                                                                                                                                                                                                                      | part of annotation dbs | https://faylward.github.io/GVDB/                               | 2021             | https://doi.org/10.1101/2021.05.05.442809                              |
| MetHMMDB               | Вatabase containing HMM profiles representing microbial metal resistance genes                                                                                                                                                                                                                                            | part of annotation dbs | https://github.com/kciuchcinski/MetHMMDB                       | 2025             | https://www.biorxiv.org/content/10.1101/2024.12.26.629440v2            |
| NMPFamsDB              | Novel Metagenome Protein Families Database                                                                                                                                                                                                                                                                                | part of annotation dbs | https://www.nmpfamsdb.org/                                     | 2024             | https://doi.org/10.1093/nar/gkad800                              |


### Download optional databases depending on choice of tools to use

#### Taxonomy assignment

| Database           | Notes                                                                                                    | Input argument                   | Version | Link                                                                                                                | Publication                                        |
| ------------------ | -------------------------------------------------------------------------------------------------------- |----------------------------------|---------|---------------------------------------------------------------------------------------------------------------------|----------------------------------------------------|
| geNomad            | geNomad database                                                                                         | `--genomad_db`                   | 4.1.1   | https://github.com/apcamargo/genomad#downloading-the-database                                                       | https://www.nature.com/articles/s41587-023-01953-y |
| ViPhOGs            | used for taxonomy assignment                                                                             | `--viphog_db`                    | 2021    | ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/viral-pipeline/hmmer*databases/vpHMM_database*{VERSION}.tar.gz       | https://www.mdpi.com/1999-4915/13/6/1164           |
| ViPhOG metadata    | metadata for filtering ViPhOGs according to taxonomy updates by the [ICTV](https://ictv.global/taxonomy) | `--additional_viphog_model_data` | 2021    | ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/viral-pipeline/additional*data_vpHMMs*{VERSION}.tsv                  | -                                                  |
| ViPhOG factor file | requred for taxonomy assignment                                                                          | `--viphogs_factor_file`          | 2021    | https://github.com/EBI-Metagenomics/emg-viral-pipeline/blob/master/references/viphogs_cds_per_taxon_cummulative.csv | -                                                  |
| ViPhOG NCBI db     | requred for taxonomy assignment                                                                          | `--ncbi_db`                      | 2023    | ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/viral-pipeline/2022-11-01_ete3_ncbi_tax.sqlite.gz                    | -                                                  |
| VITAP              | VITAP database                                                                                           | `--vitap_db`                     | 1.12    | https://github.com/DrKaiyangZheng/VITAP#use-the-pre-built-database                                                  | https://www.nature.com/articles/s41467-025-57500-7 |

#### Host detection

| Database | Notes          | Input argument | Link                                                         | Version            | Publication                                        |
| -------- | -------------- | -------------- | ------------------------------------------------------------ |--------------------|----------------------------------------------------|
| iPHoP    | iPHoP database | `--iphop_db`   | https://github.com/simroux/iphop#host-databases-and-versions | iPHoP_db_Jun25_rw  | https://pmc.ncbi.nlm.nih.gov/articles/PMC10155999/ |

#### Antimicrobial resistance detection

| Database  | Notes                                                                                                      | Input argument       | Link                                                  | Version            | Publication                                           |
| --------- | ---------------------------------------------------------------------------------------------------------- | -------------------- | ----------------------------------------------------- |--------------------|-------------------------------------------------------|
| AMRFinder | AMRFinderPlus database                                                                                     | `--amrfinderplus_db` | https://github.com/ncbi/amr/wiki/New-in-AMRFinderPlus | 4.0 (2025-07-16.1) | https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8208984/ |
| DeepARG   | DeepARG database. Also check arguments `--deeparg_db_version`, `--deeparg_model`, `--deeparg_tool_version` | `--deeparg_db`       | https://github.com/gaarangoa/deeparg                  | 2 ( tool v1.0.4)   | https://doi.org/10.1186/s40168-018-0401-z             |
| RGI       | RGI database                                                                                               | `--rgi_db`           | https://github.com/arpcard/rgi                        | 4.0.1              | https://pubmed.ncbi.nlm.nih.gov/36263822/             |

## Databases for Plasmids

#### PlaSquid

Download databases for PlaSquid workflow from https://github.com/mgimenez720/plaSquid/tree/master/data (last time updated 17 Feb 2025).

| Database file   | Argument         |
| --------------- | ---------------- |
| All_Rep.hmm     | `--repsearch_db` |
| Arq_RIP_new.RDS | `--repfilter_db` |
| All_RNA_Inc.hmm | `--rna_inc_db`   |
| All_Inc.hmm     | `--incsearch_db` |
| All_MOB.hmm     | `--mobsearch_db` |

#### MOB-suite typer

Download and prepare database with `mob_init`. Documentation: https://github.com/phac-nml/mob-suite#mob-init
