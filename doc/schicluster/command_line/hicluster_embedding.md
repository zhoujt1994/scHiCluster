# hicluster embedding

## Command Docs
```bash
usage: hicluster embedding [-h] --cell_table_path CELL_TABLE_PATH --output_dir
                           OUTPUT_DIR [--chrom_size_path CHROM_SIZE_PATH]
                           [--dim DIM] [--dist DIST] [--resolution RESOLUTION]
                           [--scale_factor SCALE_FACTOR] [--cpu CPU]
                           [--norm_sig] [--save_model] [--save_raw]

optional arguments:
  -h, --help            show this help message and exit
  --chrom_size_path CHROM_SIZE_PATH
                        Path to chromosome sizes file (default: None)
  --dim DIM
  --dist DIST
  --resolution RESOLUTION
  --scale_factor SCALE_FACTOR
  --cpu CPU
  --norm_sig
  --save_model
  --save_raw

required arguments:
  --cell_table_path CELL_TABLE_PATH
                        Path to cell table file (default: None)
  --output_dir OUTPUT_DIR
                        Path to embedding output dir (default: None)
```