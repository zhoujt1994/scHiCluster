# hicluster domain

## Command Docs
```bash
usage: hicluster domain [-h] --cell_table_path CELL_TABLE_PATH --output_prefix
                        OUTPUT_PREFIX [--resolution RESOLUTION]
                        [--window_size WINDOW_SIZE] [--cpu CPU]

optional arguments:
  -h, --help            show this help message and exit
  --resolution RESOLUTION
                        Matrix resolution (default: 25000)
  --window_size WINDOW_SIZE
                        Window size for calculating insulation score (default:
                        10)
  --cpu CPU             Number of CPUs to use (default: 10)

required arguments:
  --cell_table_path CELL_TABLE_PATH
                        Path to cell table file (default: None)
  --output_prefix OUTPUT_PREFIX
                        Output files prefix. The domain boundary matrix will
                        be saved as {output_prefix}.boundary.h5ad
                        (anndata.AnnData); The insulation score matrix will be
                        saved as {output_prefix}.insulation.nc
                        (xarray.DataSet) (default: None)
```