# Prepare and run imputation

We noticed strong signals at ENCODE blacklist regions that looks artificial. Therefore we first filter the chromatin contacts by the blacklist. In snm3C-seq data, we also provide a 2d blacklist to remove other fake signals.
```bash
# Generate a contact file list
ls raw/*/contacts_pairs/*.contacts.pairs.txt.gz | awk -v p=$PWD '{printf("%s/%s\n", p, $0)}' > contact_table.txt
paste <(awk -F'/' '{print $NF}' contact_table.txt | cut -d. -f1) contact_table.txt | sort -k1,1 > contact_table.tsv

# Remove blacklist
hicluster filter-contact --output_dir rmbkl/ --blacklist_1d_path mm10-blacklist.v2.bed.gz --chr1 1 --pos1 2 --chr2 3 --pos2 4 --contact_table contact_table.tsv --chrom_size_path chrom_sizes.txt
```

The following command will generate a snakemake_cmd.txt file, where each line is a snakemake command to run imputation with certain parameters for a number of cells (batch_size). Thus, you can append each line to your job submission template to parallelize between batches. The parallelization between cells in the same batch will be handled by snakemake on a single node. You can decide the batch_size based on your super computer setting. In the example below, our imputation would be run on nodes using 96 cpus per node, and we assigned 1536 cells to a single job so that it would not create too many jobs to avoid long queue time.
```bash
# Generate a contact file list
ls rmbkl/* | awk -v p=$PWD '{printf("%s/%s\n", p, $0)}' > contact_table_rmbkl.txt
paste <(awk -F'/' '{print $NF}' contact_table_rmbkl.txt | cut -d. -f1) contact_table_rmbkl.txt | sort -k1,1 > contact_table_rmbkl.tsv

# Prepare imputation 
# Contact file list are divided into smaller chunks assigned to different jobs
# Commands to run imputation are generated
hicluster prepare-impute --cell_table contact_table_rmbkl.tsv --batch_size 1536 --pad 1 --cpu_per_job 96 --chr1 1 --pos1 2 --chr2 3 --pos2 4 --output_dir impute/100K/ --chrom_size_path /anvil/scratch/x-zhou/Tan2021/chrom_sizes.txt --output_dist 500000000 --window_size 500000000 --step_size 500000000 --resolution 100000
hicluster prepare-impute --cell_table contact_table_rmbkl.tsv --batch_size 1536 --pad 2 --cpu_per_job 96 --chr1 1 --pos1 2 --chr2 3 --pos2 4 --output_dir impute/25K/ --chrom_size_path /anvil/scratch/x-zhou/Tan2021/chrom_sizes.txt --output_dist 5050000 --window_size 500000000 --step_size 500000000 --resolution 25000
hicluster prepare-impute --cell_table contact_table_rmbkl.tsv --batch_size 1536 --pad 2 --cpu_per_job 96 --chr1 1 --pos1 2 --chr2 3 --pos2 4 --output_dir impute/10K/ --chrom_size_path /anvil/scratch/x-zhou/Tan2021/chrom_sizes.txt --output_dist 5050000 --window_size 30000000 --step_size 10000000 --resolution 10000
mkdir -p sbatch/impute/
cat impute/*/snakemake_cmd.txt > sbatch/impute/snakemake_cmd.txt
```
Then append each line of snakemake_cmd.txt to your job submission template. Submit the jobs and run imputation! With our setting (standard node on anvil.rcac.purdue.edu), imputation of 1536 cells at 100kb resolution takes 30 minutes, 25kb resolution takes 2.5 hours, and 10kb resolution takes 5 hours.  

After Imputation, you can use the following commands to generate single-cell embedding and several cell-by-feature matrices, including compartment score matrix, insulation score matrix, domain boundary matrix, and gene interaction strength matrix. These commands are also resource consuming, so please consider submitting them to computing nodes.
```bash
# Compute embedding
ls impute/100K/chunk*/*.cool | awk -v p=$PWD '{printf("%s/%s\n", p, $0)}' > impute/100K/cell_table.txt
paste <(awk -F'/' '{print $NF}' impute/100K/cell_table.txt | cut -d. -f1) impute/100K/cell_table.txt | sort -k1,1 > impute/100K/cell_table.tsv
hicluster embedding --cell_table_path /anvil/scratch/x-zhou/Tan2021/impute/100K/cell_table.tsv --output_dir /anvil/scratch/x-zhou/Tan2021/dataset/embedding --dim 50 --dist 1000000 --resolution 100000 --scale_factor 100000 --norm_sig --save_raw --cpu 20

# Compute compartment score using raw contact matrices
## CpG density
hicluster cpg-ratio --fasta_path /anvil/projects/x-mcb130189/mm10/mm10_with_chrl.fa --hdf_output_path cpg_ratio_100k.hdf --chrom_size_path /anvil/scratch/x-zhou/Tan2021/chrom_sizes.txt --resolution 100000
## Using raw matrices
hicluster compartment --cell_table_path /anvil/scratch/x-zhou/Tan2021/contact_table_rmbkl.tsv --output_prefix /anvil/scratch/x-zhou/Tan2021/dataset/Tan2021.raw --cpg_profile_path cpg_ratio_100k.hdf --cpu 96 --resolution 100000 --chr1 1 --pos1 2 --chr2 3 --pos2 4 --chrom_size_path /anvil/scratch/x-zhou/Tan2021/chrom_sizes.txt
## Using imputed matrices
hicluster compartment --cell_table_path /anvil/scratch/x-zhou/Tan2021/impute/100K/cell_table.tsv --output_prefix /anvil/scratch/x-zhou/Tan2021/dataset/Tan2021.impute --cpg_profile_path cpg_ratio_100k.hdf --cpu 96

# Compute insulation score and domain boundary
ls impute/25K/chunk*/*.cool | awk -v p=$PWD '{printf("%s/%s\n", p, $0)}' > impute/25K/cell_table.txt
paste <(awk -F'/' '{print $NF}' impute/25K/cell_table.txt | cut -d. -f1) impute/25K/cell_table.txt | sort -k1,1 > impute/25K/cell_table.tsv
hicluster domain --cell_table_path /anvil/scratch/x-zhou/Tan2021/impute/25K/cell_table.tsv --output_prefix /anvil/scratch/x-zhou/Tan2021/dataset/Tan2021 --resolution 25000 --window_size 10 --cpu 96

# Compute gene contact score
ls impute/10K/chunk*/*.cool | awk -v p=$PWD '{printf("%s/%s\n", p, $0)}' > impute/10K/cell_table.txt
paste <(awk -F'/' '{print $NF}' impute/10K/cell_table.txt | cut -d. -f1) impute/10K/cell_table.txt | sort -k1,1 > impute/10K/cell_table.tsv
hicluster gene-score --cell_table /anvil/scratch/x-zhou/Tan2021/impute/10K/cell_table.tsv --gene_meta /anvil/projects/x-mcb130189/gene/vm22/gencode.vM22.annotation.gene.sorted.bed.gz --res 10000 --output_hdf /anvil/scratch/x-zhou/Tan2021/dataset/Tan2021.geneimputescore.hdf --chrom_size /anvil/scratch/x-zhou/Tan2021/chrom_sizes.txt --cpu 96 --mode impute
```
