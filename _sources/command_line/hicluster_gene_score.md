# hicluster gene-score
This command generate cell by gene hdf matrix.

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
                        Contain all the cool file information in twotab-
                        separated columns: 1. cell_uid, 2. file_path. No
                        header (default: None)
  --gene_meta_path GENE_META_PATH
                        Contain all gene information in four tab-seperated
                        columns: 1. chromosome, 2. start, 3. end, 4. gene_id.
                        No header (default: None)
  --resolution RESOLUTION
                        Resolution of cool file; normally use resolution at
                        10k (default: 10000)
  --output_hdf_path OUTPUT_HDF_PATH
                        Full path to output file (default: None)
  --chrom_size_path CHROM_SIZE_PATH
                        Path to UCSC chrom size file. Contain all the
                        chromosome information in two tab-separated columns:
                        1. chromosome name, 2. chromosome length. No header
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

## Command Breakdown
```bash
--cell_table_path impute/10K/cell_table.tsv
```
Specify the file paths of the imputed cool files in this line(e.g. /home/qzeng_salk_edu/project/aging/230711_m3C/impute/10K/chunk0/AMB_220712_18mo_12D_13B_2_P4-1-I15-G2.cool). Here is an example of what the contact table looks like:

```bash
cell_1  imputed_hic_cool_path_1
cell_2  imputed_hic_cool_path_2
cell_3  imputed_hic_cool_path_3
```
The first column indicates the cell name (e.g. AAMB_220712_18mo_12D_13B_2_P4-1-I15-G2) whereas the second column indicates the imputed cool file path of the cell. Make sure the two parts are separated by a tab; also make sure the file has no header.

The output file is a cell by gene matrix, values indicating contact probability on each gene.