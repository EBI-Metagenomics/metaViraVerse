# EBI-Metagenomics/metaviraverse: MGnify internal Documentation

> [!NOTE]
> That documentation page is only for EBI cluster users

## Building input dataset

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
