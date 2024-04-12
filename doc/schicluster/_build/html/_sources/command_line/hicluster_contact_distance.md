# hicluster contact-distance
This step calculate the contacts number in different genomic distances and the sparsity of contact matrices at certain resolution for all chromosomes.

## Command Docs
```bash
usage: hicluster contact-distance [-h] --contact_table CONTACT_TABLE
                                  --chrom_size_path CHROM_SIZE_PATH
                                  --output_prefix OUTPUT_PREFIX
                                  [--resolution RESOLUTION] [--chr1 CHROM1]
                                  [--pos1 POS1] [--chr2 CHROM2] [--pos2 POS2]
                                  [--cpu CPU]

optional arguments:
  -h, --help            show this help message and exit
  --resolution RESOLUTION
                        Resolution of contact length (default: 10000)
  --chr1 CHROM1         0 based index of chr1 column. (default: 1)
  --pos1 POS1           0 based index of pos1 column. (default: 2)
  --chr2 CHROM2         0 based index of chr2 column. (default: 5)
  --pos2 POS2           0 based index of pos2 column. (default: 6)
  --cpu CPU             number of cpus to parallel. (default: 20)

required arguments:
  --contact_table CONTACT_TABLE
                        Contain all the cell contact file after blacklist
                        region removwl; information in two tab-separated
                        columns: 1. cell_uid, 2. file_path. No header
                        (default: None)
  --chrom_size_path CHROM_SIZE_PATH
                        Path to UCSC chrom size fileContain all the
                        chromosomeinformation in two tab-separated columns:
                        1.chromosome name, 2. chromosome length. No header
                        (default: None)
  --output_prefix OUTPUT_PREFIX
                        Output hdf file prefix including the directory
                        (default: None)

```

## Command Example
```bash
hicluster contact-distance \
--contact_table contact_table_rmbkl.tsv \
--chrom_size_path /data/aging/ref/m3C/mm10.main.nochrM.nochrY.chrom.sizes \
--output_prefix contact_distance \
--resolution 10000 \
--chr1 1 \
--pos1 2 \
--chr2 5 \
--pos2 6 \
--cpu 20
```

## Command Break Down
```bash
--cell_table contact_table_rmbkl.tsv
```
Specify the file paths of the contact files after removing blacklist regions in this line(e.g. /home/qzeng_salk_edu/project/aging/230711_m3C/rmbkl/AMB_220712_18mo_12D_13B_2_P4-1-I15-K1.contact.rmbkl.tsv.gz). Here is an example of what the contact_table_rmbkl.tsv looks like

```bash
cell_1 absolute_hic_rmbkl_contact_path_1
cell_2 absolute_hic_rmbkl_contact_path_2
cell_3 absolute_hic_rmbkl_contact_path_3
```
The first column indicates the cell name (e.g. AMB_220712_18mo_12D_13B_2_P4-1-I15-K1) whereas the second column indicates the HiC contact file path after removing blacklist of the cell. Make sure the two parts are separated by a tab; also make sure the file has no header.

The output file of this command are contact_distance_decay.hdf5 and contact_distance_chromsparsity.hdf5, which can be read using pd.read_hdf coomand. The decay file records the number of contacts in different genomic distances, while the chromsparsity file shows the total number of contacts on each chromosome. 