# hicluster contact-distance

## Command Docs
```bash
usage: hicluster contact-distance [-h] [--contact_table CONTACT_TABLE]
                                  [--chrom_size_path CHROM_SIZE_PATH]
                                  [--output_prefix OUTPUT_PREFIX]
                                  [--resolution RESOLUTION] [--chr1 CHROM1]
                                  [--pos1 POS1] [--chr2 CHROM2] [--pos2 POS2]
                                  [--cpu CPU]

optional arguments:
  -h, --help            show this help message and exit
  --resolution RESOLUTION
                        Resolution of cool file (default: 10000)
  --chr1 CHROM1         0 based index of chr1 column. (default: 1)
  --pos1 POS1           0 based index of pos1 column. (default: 2)
  --chr2 CHROM2         0 based index of chr2 column. (default: 5)
  --pos2 POS2           0 based index of pos2 column. (default: 6)
  --cpu CPU             number of cpus to parallel. (default: 20)

required arguments:
  --contact_table CONTACT_TABLE
                        Full path to cell contact files (default: None)
  --chrom_size_path CHROM_SIZE_PATH
                        Chromsome size file with only chromosomes to use
                        (default: None)
  --output_prefix OUTPUT_PREFIX
                        Output hdf file prefix including the directory
                        (default: None)
```