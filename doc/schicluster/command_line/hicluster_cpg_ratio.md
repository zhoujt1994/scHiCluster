# hicluster cpg-ratio
This step willl calculate the CpG ratio of each chromotin regions at specific resolution. 

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
                        Path to UCSC chrom size fileContain all the
                        chromosomeinformation in two tab-separated columns:1.
                        chromosome name, 2. chromosome length. No header
                        (default: None)
  --dim DIM
  --dist DIST
  --resolution RESOLUTION
                        Resolution for embedding. Consistent with resolution
                        of imputed contact files (default: 100000)
  --scale_factor SCALE_FACTOR
  --cpu CPU
  --norm_sig
  --save_model
  --save_raw

required arguments:
  --cell_table_path CELL_TABLE_PATH
                        Contain all the imputed contact files information in
                        twotab-separated columns: 1. cell_uid, 2. file_path.
                        No header (default: None)
  --output_dir OUTPUT_DIR
                        Path to the output directory of the embedding output
                        (default: None)
(allcools) ➜  scHiCluster git:(rachel-hicluster) ✗ hicluster cpg-ratio -h
usage: hicluster cpg-ratio [-h] --fasta_path FASTA_PATH --hdf_output_path
                           HDF_OUTPUT_PATH [--cell_url CELL_URL]
                           [--chrom_size_path CHROM_SIZE_PATH]
                           [--resolution RESOLUTION]

optional arguments:
  -h, --help            show this help message and exit
  --cell_url CELL_URL   Path to a cell Cooler URL (default: None)
  --chrom_size_path CHROM_SIZE_PATH
                        Path to UCSC chrom size fileContain all the
                        chromosomeinformation in two tab-separated columns:1.
                        chromosome name, 2. chromosome length. No header
                        (default: None)
  --resolution RESOLUTION

required arguments:
  --fasta_path FASTA_PATH
                        Path to genome FASTA file (default: None)
  --hdf_output_path HDF_OUTPUT_PATH
                        Output path of the CpG ratio (default: None)
```

## Command Example
Here is an example of calculating CpG ratio at each Chrom100K (resolution 100000). 

```bash
hicluster cpg-ratio \
--fasta_path /data/aging/ref/m3C/mm10_with_chrl.fa \
--hdf_output_path cpg_ratio_100k.hdf \
--chrom_size_path /data/aging/ref/m3C/mm10.main.nochrM.nochrY.chrom.sizes \
--resolution 100000
```
The code belowe will generate outfile named cpg_ratio_100k.hdf, which will be needed in computing compartment scroe as a way tp normalized CpG backgroud. 