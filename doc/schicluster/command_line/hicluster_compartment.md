# hicluster compartment
This step shows how we compute compartment score. 

## Command Docs
```bash
usage: hicluster compartment [-h] --cell_table_path CELL_TABLE_PATH
                             --output_prefix OUTPUT_PREFIX --cpg_profile_path
                             CPG_PROFILE_PATH [--cpu CPU] [--calc_strength]
                             [--mode {tsv,cool}]
                             [--chrom_size_path CHROM_SIZE_PATH]
                             [--resolution RESOLUTION] [--chr1 CHROM1]
                             [--chr2 CHROM2] [--pos1 POS1] [--pos2 POS2]

optional arguments:
  -h, --help            show this help message and exit
  --cpu CPU             Number of CPUs to use (default: 10)
  --calc_strength       Calculate compartment strength summary (default:
                        False)
  --mode {tsv,cool}     cool or tsv (default: cool)
  --chrom_size_path CHROM_SIZE_PATH
                        Path to UCSC chrom size fileContain all the chromosome
                        information in two tab-separated columns: 1.
                        chromosome name, 2. chromosome length. No header
                        (default: None)
  --resolution RESOLUTION
  --chr1 CHROM1         0 based index of chr1 column. (default: 1)
  --chr2 CHROM2         0 based index of chr2 column. (default: 5)
  --pos1 POS1           0 based index of pos1 column. (default: 2)
  --pos2 POS2           0 based index of pos2 column. (default: 6)

required arguments:
  --cell_table_path CELL_TABLE_PATH
                        Contain all the cell contact file information in two
                        tab-separated columns: 1. cell_uid, 2. file_path. No
                        header (default: None)
  --output_prefix OUTPUT_PREFIX
                        Output files prefix. The compartment score matrix will
                        be saved as {output_prefix}.compartment.nc (use
                        xr.open_dataset() to load). (default: None)
  --cpg_profile_path CPG_PROFILE_PATH
                        Genome bins CpG ratio. Use "schicluster cpg-ratio" to
                        calculate (default: None)

```

## Command Example
This will compute compartment score using raw contact files.
```bash
hicluster compartment \
--cell_table_path contact_table_rmbkl.tsv \
--output_prefix  dataset/raw \
--cpg_profile_path cpg_ratio_100k.hdf \
--cpu 48 \
--resolution 100000 \
--chr1 1 \
--pos1 2 \
--chr2 5 \
--pos2 6 \
--chrom_size_path /data/aging/ref/m3C/mm10.main.nochrM.nochrY.chrom.sizes \
--mode tsv
```

This will compute compartment score using imputed contact files.
```bash
hicluster compartment \
--cell_table_path impute/100K/cell_table.tsv \
--output_prefix dataset/impute \
--cpg_profile_path cpg_ratio_100k.hdf \
--chr1 1 \
--pos1 2 \
--chr2 5 \
--pos2 6 \
--cpu 48
```

## Command Break Down
```bash
--cpg_profile_path cpg_ratio_100k.hdf
```
The cpg_ratio_100k.hdf will be generated by hicluster cpg_ratio command.

```bash
--cell_table_path contact_table_rmbkl.tsv
--cell_table_path impute/100K/cell_table.tsv
```
The first column indicates the cell name (e.g. AMB_220712_18mo_12D_13B_2_P4-1-I15-A13) whereas the second column indicates the path to hic contact file after blacklist removal or path to cool cile. Make sure the two parts are separated by a tab; also make sure the file has no header and index.

The output file is an xarray file which can be read load using xrray.open_dataset().