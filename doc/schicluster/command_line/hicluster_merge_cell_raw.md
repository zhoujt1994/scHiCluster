# hicluster merge-cell-raw
This step merge single-cell contacts by summing up. 

## Command Docs
```bash
usage: hicluster merge-cell-raw [-h] --cell_table CELL_TABLE --chrom_size_path
                                CHROM_SIZE_PATH --output_file OUTPUT_FILE
                                [--resolution RESOLUTION] [--chr1 CHROM1]
                                [--pos1 POS1] [--chr2 CHROM2] [--pos2 POS2]
                                [--min_pos_dist MIN_POS_DIST]

optional arguments:
  -h, --help            show this help message and exit
  --resolution RESOLUTION
                        Resolution of cool file (default: 5000)
  --chr1 CHROM1         0 based index of chr1 column. (default: 1)
  --pos1 POS1           0 based index of pos1 column. (default: 2)
  --chr2 CHROM2         0 based index of chr2 column. (default: 5)
  --pos2 POS2           0 based index of pos2 column. (default: 6)
  --min_pos_dist MIN_POS_DIST
                        Minimum distance for a fragment to be considered.
                        (default: 2500)

required arguments:
  --cell_table CELL_TABLE
                        Contain all the cell contact file after
                        blacklistremoval in two tab-separated columns: 1.
                        cell_uid, 2.file_path. No header (default: None)
  --chrom_size_path CHROM_SIZE_PATH
                        Path to UCSC chrom size fileContain all the
                        chromosomeinformation in two tab-separated columns: 1.
                        chromosome name, 2. chromosome length. No header
                        (default: None)
  --output_file OUTPUT_FILE
```

## Command Examples
```bash
hicluster merge-cell-raw \
--cell_table contact_table_rmbkl.tsv \
--chrom_size_path /data/aging/ref/m3C/mm10.main.nochrM.nochrY.chrom.sizes \
--output_file dataset/merge-cell-raw.5kb \
--resolution 5000 \
--chr1 1 \
--pos1 2 \
--chr2 5 \
--pos2 6 \
--min_pos_dist 2500
```