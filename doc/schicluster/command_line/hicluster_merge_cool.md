# hicluster merge-cool
This step merge single-cell cool files or contact files by averaging over cells. s

## Command Docs
```bash
usage: hicluster merge-cool [-h] --input_cool_tsv_file INPUT_COOL_TSV_FILE
                            --output_cool OUTPUT_COOL

optional arguments:
  -h, --help            show this help message and exit

required arguments:
  --input_cool_tsv_file INPUT_COOL_TSV_FILE
                        Contain all the imputed cool file information in a tsv
                        file: 1. file_path. No header (default: None)
  --output_cool OUTPUT_COOL
                        Full path to output merged cool file (default: None)
```

## Command Examples
```bash
hicluster merge-cool \
--input_cool_tsv_file impute/100K/cell_table_2.tsv \
--output_cool dataset/merge-cool.100kb
```

