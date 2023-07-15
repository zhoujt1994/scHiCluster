# hicluster merge-cell-raw

## Command Docs
```bash
usage: hicluster merge-cell-raw [-h] [--cell_table CELL_TABLE]
                                [--chrom_size_path CHROM_SIZE_PATH]
                                [--output_file OUTPUT_FILE]
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
                        Full path to cell cool table (default: None)
  --chrom_size_path CHROM_SIZE_PATH
                        Chromsome size file with only chromosomes to use
                        (default: None)
  --output_file OUTPUT_FILE
                        Full path to output file (default: None)
```