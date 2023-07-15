# hicluster prepare-impute

After removing the blacklist regions, the following command will generate a file called "snakemake_cmd.txt". Each line in this file represents a Snakemake command for running imputation with specific parameters for a certain number of cells (batch_size). By appending each line to your job submission template, you can parallelize the imputation process between batches. The parallelization within each batch, handling multiple cells, will be managed by Snakemake on a single node.

The choice of batch_size should be based on the settings of your supercomputer. Previsou benchmark indicate that, jobs executed on nodes with 96 CPUs per node. We assign 1536 cells to each job in order to strike a balance and avoid creating an excessive number of jobs, which could result in long queue times.

## Command Docs
```bash
usage: hicluster prepare-impute [-h] --output_dir OUTPUT_DIR --chrom_size_path
                                CHROM_SIZE_PATH --output_dist OUTPUT_DIST
                                --window_size WINDOW_SIZE --step_size
                                STEP_SIZE --resolution RESOLUTION
                                [--input_scool INPUT_SCOOL]
                                [--cell_table CELL_TABLE]
                                [--batch_size BATCH_SIZE] [--logscale]
                                [--pad PAD] [--std STD] [--rp RP] [--tol TOL]
                                [--min_cutoff MIN_CUTOFF]
                                [--cpu_per_job CPU_PER_JOB] [--chr1 CHROM1]
                                [--chr2 CHROM2] [--pos1 POS1] [--pos2 POS2]

optional arguments:
  -h, --help            show this help message and exit
  --input_scool INPUT_SCOOL
                        Path to input scool file (default: None)
  --cell_table CELL_TABLE
                        Contain all the cell contact file after blacklist
                        removalin two tab-separated columns: 1. cell_uid, 2.
                        file_path.No header (default: None)
  --batch_size BATCH_SIZE
                        Number of cells to include in each batch run (default:
                        100)
  --logscale
  --pad PAD
  --std STD
  --rp RP
  --tol TOL
  --min_cutoff MIN_CUTOFF
  --cpu_per_job CPU_PER_JOB
                        Number of cpus to parallel. (default: 10)
  --chr1 CHROM1         0 based index of chr1 column. (default: 1)
  --chr2 CHROM2         0 based index of chr2 column. (default: 5)
  --pos1 POS1           0 based index of pos1 column. (default: 2)
  --pos2 POS2           0 based index of pos2 column. (default: 6)

required arguments:
  --output_dir OUTPUT_DIR
                        Path to output directory containing the single-cell
                        cool files after imputation (default: None)
  --chrom_size_path CHROM_SIZE_PATH
                        Path to UCSC chrom size fileContain all the chromosome
                        information in two tab-separated columns: 1.
                        chromosome name, 2. chromosome length. No header
                        (default: None)
  --output_dist OUTPUT_DIST
                        Maximum distance for a contact to be written to the
                        output (default: 500000000)
  --window_size WINDOW_SIZE
                        The size in base pairs of sliding window for
                        imputation (default: 500000000)
  --step_size STEP_SIZE
                        The step size in base pairs that the sliding window
                        moves (default: 10000000)
  --resolution RESOLUTION
                        Resolution for imputation (default: None)
```

## Eaxmple 1: Imputation at 100kb resolution 
```bash 
hicluster prepare-impute 
    --cell_table contact_table_rmbkl.tsv \
    --batch_size 1536 \
    --pad 1 \
    --cpu_per_job 30 \
    --chr1 1 \
    --pos1 2 \
    --chr2 5 \
    --pos2 6 \
    --output_dir impute/100K/ \
    --chrom_size_path /data/aging/ref/m3C/mm10.main.nochrM.nochrY.chrom.sizes \
    --output_dist 500000000 \
    --window_size 500000000 \
    --step_size 500000000 \
    --resolution 100000 \
```

```bash
--cell_table contact_table_rmbkl.tsv
```
Specify the file paths of the contact files after removing blacklist regions in this line(e.g. /home/qzeng_salk_edu/project/aging/230711_m3C/rmbkl/AMB_220712_18mo_12D_13B_2_P4-1-I15-K1.contact.rmbkl.tsv.gz). Here is an example of what the contact_table_rmbkl.tsv looks like

```bash
cell_1 absolute_hic_rmbkl_contact_path_1
cell_2 absolute_hic_rmbkl_contact_path_2
cell_3 absolute_hic_rmbkl_contact_path_3
```
The first column indicates the cell name (e.g. AMB_220712_18mo_12D_13B_2_P4-1-I15-K1) whereas the second column indicates the HiC contact file path after removing blacklist of the cell. Make sure the two parts are separated by a tab; also make sure the file has no header and index.


You don't need to create the folder before running.
```bash
--output_dir impute/100K/
```

After running the above command, we get a file names snakemake_cmd.txt. Run the command to get imputed cool files at 100kb resolution :
```bash
bash /home/qzeng_salk_edu/project/aging/230711_m3C/impute/100K/snakemake_cmd.txt 
```


## Eaxmple 2: Imputation at 25kb resolution
```bash
hicluster prepare-impute 
    --cell_table contact_table_rmbkl.tsv \
    --batch_size 1536 \
    --pad 1 \
    --cpu_per_job 30 \
    --chr1 1 \
    --pos1 2 \
    --chr2 5 \
    --pos2 6 \
    --output_dir impute/100K/ \
    --chrom_size_path /data/aging/ref/m3C/mm10.main.nochrM.nochrY.chrom.sizes \
    --output_dist 5050000 \
    --window_size 500000000 \
    --step_size 500000000 \
    --resolution 25000 \
```

## Eaxmple 3: Imputation at 10kb resolution
```bash
hicluster prepare-impute 
    --cell_table contact_table_rmbkl.tsv \
    --batch_size 1536 \
    --pad 1 \
    --cpu_per_job 30 \
    --chr1 1 \
    --pos1 2 \
    --chr2 5 \
    --pos2 6 \
    --output_dir impute/100K/ \
    --chrom_size_path /data/aging/ref/m3C/mm10.main.nochrM.nochrY.chrom.sizes \
    --output_dist 5050000 \
    --window_size 30000000 \
    --step_size 10000000 \
    --resolution 10000 \
```




