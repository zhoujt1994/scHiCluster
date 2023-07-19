# hicluster embedding
This step will generate cell by 100kb-pair contact matrix which can be used for cell embedding

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
                        Resolution for embedding.Consistent with resolution of
                        imputed contact files (default: 100000)
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

```

## Command Example
Here is an example to generate chrom by chrom contact matrix. 
```
hicluster embedding \
    --cell_table_path cell_table.tsv \
    --output_dir dataset/embedding \
    --dim 50 \
    --dist 1000000 \
    --resolution 100000 \
    --scale_factor 100000 \
    --norm_sig \
    --save_raw \
    --cpu 20 
```

## Command Breakdown

```bash
--cell_table_path cell_table.tsv
```
Specify the file paths of the cool files after imputtaion in this line(e.g. /home/qzeng_salk_edu/project/aging/230711_m3C/impute/100K/chunk0/AMB_220712_18mo_12D_13B_2_P4-1-I15-G2.cool). Here is an example of what the contact table looks like:

```bash
cell_1 absolute_cool_file_path_1
cell_2 absolute_cool_file_path_2
cell_3 absolute_cool_file_path_3
```
The first column indicates the cell name (e.g. AMB_220712_18mo_12D_13B_2_P4-1-I15-G2) whereas the second column indicates the cool file path of the cell. Make sure the two parts are separated by a tab; also make sure the file has no header.

```bash
--output_dir dataset/embedding
```
This is the path to the output folder for your output files and you don't need to create the folder before running. This command will save output files in output_dir/raw and output_dir/decomp. 

The output_dir/raw: npz files of each chromosome, which contains the information of cell x 100kb-pair contacts matrix (e.g. chr1.npz). 

The folder "output_dir/decomp": the concatenated contacts of all chromosomes after performing singular value decomposition (SVD) on each chromosome. (total_chrom_decomp_concat.npz). Additionally, the concatenated decomposition matrices of all chromosomes are further subjected to another round of SVD. (total_decomp.npz)

For information regading loading npz files see [here] https://numpy.org/doc/stable/reference/generated/numpy.savez.html




