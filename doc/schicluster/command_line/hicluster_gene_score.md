# hicluster gene-score
This command compute gene contact score

## Command Docs
```bash
usage: hicluster gene-score [-h] --cell_table_path CELL_TABLE_PATH
                            --gene_meta_path GENE_META_PATH --resolution
                            RESOLUTION --output_hdf_path OUTPUT_HDF_PATH
                            --chrom_size_path CHROM_SIZE_PATH [--cpu CPU]
                            [--slop SLOP] [--mode MODE] [--chr1 CHROM1]
                            [--chr2 CHROM2] [--pos1 POS1] [--pos2 POS2]

optional arguments:
  -h, --help            show this help message and exit
  --cpu CPU             CPUs to use (default: 10)
  --slop SLOP           gene slop distance on both sides (default: 0)
  --mode MODE           raw or impute (default: impute)
  --chr1 CHROM1         0 based index of chr1 column. (default: 1)
  --chr2 CHROM2         0 based index of chr2 column. (default: 5)
  --pos1 POS1           0 based index of pos1 column. (default: 2)
  --pos2 POS2           0 based index of pos2 column. (default: 6)

required arguments:
  --cell_table_path CELL_TABLE_PATH
                        Full path to cell cool table (default: None)
  --gene_meta_path GENE_META_PATH
                        Full path to bed file with region id (default: None)
  --resolution RESOLUTION
                        Resolution of cool file (default: 10000)
  --output_hdf_path OUTPUT_HDF_PATH
                        Full path to output file (default: None)
  --chrom_size_path CHROM_SIZE_PATH
                        Chromsome size file with only chromosomes to use
                        (default: None)
```

## Command Examples
```bash
hicluster gene-score \
--cell_table_path impute/10K/cell_table.tsv \
--gene_meta_path /data/aging/ref/m3C/gencode.vM22.annotation.gene.sorted.bed.gz \
--resolution 10000 \
--output_hdf_path  geneimputescore.hdf \
--chrom_size_path /data/aging/ref/m3C/mm10.main.nochrM.nochrY.chrom.sizes \
--cpu 48 
--mode impute
```