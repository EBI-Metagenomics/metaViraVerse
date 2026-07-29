## Databases for Viruses and Prophages

### Download mandatory databases

| Database               | Notes                                       | Input argument | Link | Publication |
|------------------------|---------------------------------------------|---------------|------------------------|-------------|
| CheckV                 | CheckV database                             | `--checkv_db` | https://portal.nersc.gov/CheckV/checkv-db-v1.5.tar.gz | https://www.nature.com/articles/s41587-020-00774-7 |
| Protein annotation DBs | Databases were taken from MetaCerberus tool | `--annotation_db` | https://github.com/raw-lab/MetaCerberus#metacerberus-databases | https://academic.oup.com/bioinformatics/article/40/3/btae119/7616988 |

### Download optional databases depending on choice of tools to use

#### Taxonomy assignment
| Database           | Notes                                                                                                    | Input argument | Link | Publication |
|--------------------|----------------------------------------------------------------------------------------------------------|---------------|------------------------|-------------|
| geNomad            | geNomad database                                                                                         | `--genomad_db` | https://github.com/apcamargo/genomad#downloading-the-database | https://www.nature.com/articles/s41587-023-01953-y |
| ViPhOGs            | used for taxonomy assignment                                                                             | `--viphog_db` | ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/viral-pipeline/hmmer_databases/vpHMM_database_{VERSION}.tar.gz | https://www.mdpi.com/1999-4915/13/6/1164 |
| ViPhOG metadata    | metadata for filtering ViPhOGs according to taxonomy updates by the [ICTV](https://ictv.global/taxonomy) | `--additional_model_data` | ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/viral-pipeline/additional_data_vpHMMs_{VERSION}.tsv | - |
| ViPhOG factor file | requred for taxonomy assignment                                                               | `--factor_file` | https://github.com/EBI-Metagenomics/emg-viral-pipeline/blob/master/references/viphogs_cds_per_taxon_cummulative.csv | - |
| ViPhOG NCBI db     | requred for taxonomy assignment                                                               | `--ncbi_db` | ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/viral-pipeline/2022-11-01_ete3_ncbi_tax.sqlite.gz | - |
| VITAP              | VITAP database                                                                                | `--vitap_db` | https://github.com/DrKaiyangZheng/VITAP#use-the-pre-built-database | https://www.nature.com/articles/s41467-025-57500-7 |

#### Host detection
| Database | Notes          | Input argument | Link | Publication                                         |
|----------|----------------|---------------|------------------------|-----------------------------------------------------|
| iPHoP    | iPHoP database | `--iphop_db` | https://github.com/simroux/iphop#host-databases-and-versions| https://pmc.ncbi.nlm.nih.gov/articles/PMC10155999/  |


#### Antimicrobial resistance detection

| Database  | Notes                                                                                                      | Input argument | Link                                                   | Publication                                         |
|-----------|------------------------------------------------------------------------------------------------------------|---------------|--------------------------------------------------------|-----------------------------------------------------|
| AMRFinder | AMRFinderPlus database                                                                                     | `--amrfinderplus_db` | https://github.com/ncbi/amr/wiki/New-in-AMRFinderPlus  | https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8208984/ |
| DeepARG   | DeepARG database. Also check arguments `--deeparg_db_version`, `--deeparg_model`, `--deeparg_tool_version` | `--deeparg_db` | https://github.com/gaarangoa/deeparg  | https://doi.org/10.1186/s40168-018-0401-z |
| RGI       | RGI database                                                                                               | `--rgi_db` | https://github.com/arpcard/rgi  | https://pubmed.ncbi.nlm.nih.gov/36263822/ |

## Databases for Plasmids

#### PlaSquid
Download databases for PlaSquid workflow from https://github.com/mgimenez720/plaSquid/tree/master/data

| Database file | Argument      | 
|---------------|---------------|
| All_Rep.hmm     | `--repsearch_db` | 
| Arq_RIP_new.RDS   | `--repfilter_db` | 
| All_RNA_Inc.hmm           | `--rna_inc_db`| 
| All_Inc.hmm           | `--incsearch_db` | 
| All_MOB.hmm           | `--mobsearch_db` | 

#### MOB-suite typer

Download and prepare database with `mob_init`. Documentation: https://github.com/phac-nml/mob-suite#mob-init