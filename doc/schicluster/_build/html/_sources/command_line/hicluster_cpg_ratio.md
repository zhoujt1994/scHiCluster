# hicluster cpg-ratio
This step willl calculate the genome bins CpG ratio.

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
                        Path to UCSC chrom size fileContain all the chromosome
                        information in two tab-separated columns: 1.
                        chromosome name, 2. chromosome length. No header
                        (default: None)
  --resolution RESOLUTION
                        Resolution of the bin size in output CpG file
                        (default: 100000)

required arguments:
  --fasta_path FASTA_PATH
                        Path to genome FASTA file (default: None)
  --hdf_output_path HDF_OUTPUT_PATH
                        Output path of the CpG ratio hdf (default: None)
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
The code belowe will generate outfile named cpg_ratio_100k.hdf, which contains information of CpG ratio at each chromosome bin. This file will be needed in computing compartment score. 