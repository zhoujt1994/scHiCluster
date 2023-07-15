# hicluster cpg-ratio

## Command Docs
```bash
usage: hicluster cpg-ratio [-h] --fasta_path FASTA_PATH --hdf_output_path
                           HDF_OUTPUT_PATH [--cell_url CELL_URL]
                           [--chrom_size_path CHROM_SIZE_PATH]
                           [--resolution RESOLUTION]

optional arguments:
  -h, --help            show this help message and exit
  --cell_url CELL_URL   Path to a cell Cooler URL (default: None)
  --chrom_size_path CHROM_SIZE_PATH
                        Path to chromosome sizes file (default: None)
  --resolution RESOLUTION

required arguments:
  --fasta_path FASTA_PATH
                        Path to genome FASTA file (default: None)
  --hdf_output_path HDF_OUTPUT_PATH
                        Output path of the CpG ratio (default: None)
```