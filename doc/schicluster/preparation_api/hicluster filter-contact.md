# hicluster filter-contact
This step shows the first step for processing the HiC contact files. Here is an example of processing snm3C-seq data of mouse brain.

## Command
```bash
hicluster filter-contact \
    --output_dir rmbkl \
    --blacklist_1d_path /data/aging/ref/m3C/mm10-blacklist.v2.bed.gz \
    --blacklist_2d_path /data/aging/ref/m3C/mm10_2d_blacklist.bedpe \
    --cpu 20 \
    --chr1 1 \
    --pos1 2 \
    --chr2 5 \
    --pos2 6 \
    --contact_table contact_table.tsv \
    --chrom_size_path /data/aging/ref/m3C/mm10.main.nochrM.nochrY.chrom.sizes 
```

## Command  Breakdown

```bash
--contact_table contact_table.tsv
```
Specify the file paths of the contact files in this line(e.g. /data/AMB-F-mapping/pool_amb64/mapping_000024/mapping_000024/hic/AMB_220712_18mo_12D_13B_2_P4-1-I15-A13.hisat3n_dna.all_reads.3C.contact.tsv.gz). Here is an example of what the contact table looks like:

```bash
cell_1 absolute_hic_contact_path_1
cell_2 absolute_hic_contact_path_2
cell_3 absolute_hic_contact_path_3
```
The first column indicates the cell name (e.g. AMB_220712_18mo_12D_13B_2_P4-1-I15-A13) whereas the second column indicates the hic contact file path of the cell. Make sure the two parts are separated by a tab; also make sure the file has no header and index.

```bash
--output_dir rmbkl
```
Make sure the output directory are created before running the command. 

## Command Doc
```bash
usage: hicluster filter-contact [-h] [--contact_table CONTACT_TABLE]
                                [--chrom_size_path CHROM_SIZE_PATH]
                                [--output_dir OUTPUT_DIR]
                                [--blacklist_1d_path BLACKLIST_1D_PATH]
                                [--blacklist_2d_path BLACKLIST_2D_PATH]
                                [--blacklist_resolution RESOLUTION_2D]
                                [--not_remove_duplicates] [--chr1 CHROM1]
                                [--pos1 POS1] [--chr2 CHROM2] [--pos2 POS2]
                                [--min_pos_dist MIN_POS_DIST] [--cpu CPU]

optional arguments:
  -h, --help            show this help message and exit
  --output_dir OUTPUT_DIR
                        Output path of the contact files after blacklist
                        filtering (default: None) (default: None)
  --blacklist_1d_path BLACKLIST_1D_PATH
                        Path to blacklist region BED file, such as ENCODE
                        blacklist. Either side of the contact overlapping with
                        a blacklist region will be removed. (default: None)
  --blacklist_2d_path BLACKLIST_2D_PATH
                        Path to blacklist region pair BEDPE file. Both side of
                        the contact overlapping with the same blacklist region
                        pair will be removed. (default: None)
  --blacklist_resolution RESOLUTION_2D
                        Resolution in bps when consider the 2D blacklist
                        region pairs. (default: 10000)
  --not_remove_duplicates
                        If set, will NOT remove duplicated contacts based on
                        [chr1, pos1, chr2, pos2] values (default: True)
  --chr1 CHROM1         0 based index of chr1 column. (default: 1)
  --pos1 POS1           0 based index of pos1 column. (default: 2)
  --chr2 CHROM2         0 based index of chr2 column. (default: 5)
  --pos2 POS2           0 based index of pos2 column. (default: 6)
  --min_pos_dist MIN_POS_DIST
                        Minimum distance for a fragment to be considered.
                        (default: 0)
  --cpu CPU             number of cpus to parallel. (default: 20)

required arguments:
  --contact_table CONTACT_TABLE
                        Contain all the cell contact file information in two
                        tab-separated columns: 1. file_uid, 2. file_path. No
                        header; No index (default: None)
  --chrom_size_path CHROM_SIZE_PATH
                        Path to UCSC chrom size file. This can be generated
                        from the genome fasta or downloaded via UCSC
                        fetchChromSizes tools. (default: None)
```