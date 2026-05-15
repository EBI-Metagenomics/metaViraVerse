# EBI-Metagenomics/metaviraverse: MGnify internal Documentation

> [!NOTE]
> That documentation page is only for EBI cluster users

## Building input dataset

### MGnify genomes

#### Collect data from existing catalogues
You need to choose what catalogues you want to use and find their locations on `/nfs/public/`. \
Run fetching script [`collect_data_from_catalogues.py`](../scripts/collect_data_from_catalogues.py) (make sure you are in correct queue to access NFS)

```commandline
usage: collect_data_from_catalogues.py [-h] -p CATALOGUE_PATH [CATALOGUE_PATH ...] -o OUTPUT_PATH

Script searches for viral records in catalogue(s) GFFs and greps corresponding nucleotide and protein sequences.

options:
  -h, --help            show this help message and exit
  -p, --catalogue-path CATALOGUE_PATH [CATALOGUE_PATH ...]
                        Path to NFS location of catalogue(s)
  -o, --output-path OUTPUT_PATH
                        Path to save results (filtered gff, fna, faa)

        Script Takes as input path(s) to catalogue(s) and creates 3 files with all found viral_sequences and plasmids:
        output: catalogue_name_version.fna and catalogue_name_version.faa
```

example,

```
python3 collect_data_from_catalogues.py \
  -p nfs/catalogue_1/v1.0 nfs/catalogue_2/v1.0 \
  -o results
```

#### Outputs description:
- `<CATALOGUE_NAME>_<VERSION>_viral.fna` - FASTA file with viral sequences (viral_sequences/propahges/plasmids) found in CATALOGUE_NAME__VERSION. \
    Header format example:
    ```
    >MGYG000517142_2 prophage|311953:328572
    ```
    Where \
    `MGYG000517142` is MGnify genome identifier \
    `MGYG000517142_2` - contig identifier \
    `prophage` - sequence type \
    `311953:328572` - region coordinates \

- `<CATALOGUE_NAME>_<VERSION>_viral.faa` - FASTA file with proteins detected within chosen viral regions. \
    Header format example:
    ```
    >MGYG000517142_00528 Alpha-xylosidase
    ```
    Where \
    `MGYG000517142_00528` - protein identifier
    `Alpha-xylosidase` - protein name
- `<CATALOGUE_NAME>_<VERSION>_viral.gff` - GFF with all viral records and CDS detected within those regions \
   ```
    MGYG000517142_2 geNomad prophage        311953  328572  .       .       .       ID=MGYG000517142_2|prophage-311953:328572;mobile_element_type=prophage;taxonomy=Viruses%3BDuplodnaviria%3BHeunggongvirae%3BUroviricota%3BCaudoviricetes%3B%3B
    MGYG000517142_2 Prodigal:002006 CDS     18      1355    .       +       0       ID=MGYG000517142_00528;eC_number=3.2.1.177;Name=yicI;Dbxref=COG:COG1501;gene=yicI;inference=ab initio prediction:Prodigal:002006,similar to AA sequence:UniProtKB:P31434;locus_tag=MGYG000517142_00528;product=Alpha-xylosidase;eggNOG=1194165.CAJF01000023_gene3185;cog=G;kegg=ko:K01811;pfam=PF01055,PF21365;interpro=IPR000322,IPR013780,IPR017853,IPR048395,IPR050985
    MGYG000517142_2 Prodigal:002006 CDS     1352    2776    .       +       0       ID=MGYG000517142_00529;eC_number=3.2.1.21;inference=ab initio prediction:Prodigal:002006,similar to AA sequence:UniProtKB:P94248;locus_tag=MGYG000517142_00529;product=Bifunctional beta-D-glucosidase/beta-D-fucosidase;eggNOG=1194165.CAJF01000023_gene3186;cog=G;kegg=ko:K05350;pfam=PF00232;interpro=IPR001360,IPR017736,IPR017853
   ```
  The structure follows format: viral region and then predicted CDS inside that region. 
  
    **Important!** 

    **Region IDs** correspond to whole sequence headers in `<CATALOGUE_NAME>_<VERSION>_viral.fna`

    **CDS IDs** correspond to protein identifies from `<CATALOGUE_NAME>_<VERSION>_viral.faa`


