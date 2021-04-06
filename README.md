# scHiCluster
## Introduction
scHiCluster is a comprehensive python package for single-cell chromosome contact data analysis. It includes the identification of cell types (clusters), loop calling in cell types, and domain and compartment calling in single cells.

<img src="example/plot/Introduction.png" width="700" height="200" />  

## Installation
Running scHiCluster requries numpy, scipy, pandas, h5py, scikit-learn, opencv-python, statsmodels.  
In order to visualize the results, we also need matplotlib, umap-learn, multicore-tsne, and harmonypy to account for potential batch effects.  
First, creat a new conda environment and activate it by
```
conda create --n schicluster python==3.6.8
conda activate schicluster
```
Then install install scHiCluster by
```
pip install schicluster
```

## Usage
### General file
HiCluster requires a the chromosome size file in several following steps, where the first column is the chromosome name and the second column is the size of the chromosome in bps. The files for hg38, hg19, and mm10 are provided in file/ folder of this repo.

### Clustering
HiCluster uses linear convolution and random walk with restart to impute the chromatin contact matrices for each cell and each chromosome separately. The imputed matrices are then concatenated and used for embedding, visualation and clustering.

The input file format for scHiCluster is the sparse format contact matrices. For each cell and each chromosome, the input file should contain three columns separated by tab, representing the interacting bins and the number of reads supporting the interaction. The name of the file need to be in the format of '{cell_id}_{chromosome}.txt'.  

As an example, at 1mb resolution, if there are 10 reads in cell_1 supporting the interaction between chr1:1000000-2000000 and chr1:5000000-6000000, this should be represented as
> 1	5	10

in a single line of the file named as cell_1_chr1.txt. 

Alternatively, if you have the chromatin contact file of a single cell, the following command can be used to generate the input matrices at a specific resolution. The contact file should have the chromosome names and positions of the two interaction anchors of each contact in a single line. The option --chr1 --pos1 --chr2 --pos2 are used to indicate which columns of the file contain these information. Note that the columns should be counted start from 0. As an example, if using the [juicer-pre short format](https://github.com/aidenlab/juicer/wiki/Pre#short-format), the conmmand should be
```
hicluster generatematrix-cell --infile {input_dir}{contact_file} --outdir {output_dir} --chrom_file {chromosome_size_file} --res {resolution} --cell {cell_id} --chr1 1 --pos1 2 --chr2 5 --pos2 6
```
Then you can impute a single-cell contact matrix by
```
hicluster impute-cell --indir {raw_dir} --outdir {impute_dir}/ --cell ${cell_id} --chrom ${chromosome} --res ${resolution} --chrom_file {chromosome_size_file}
```
This can be easily parrelized across cells and chromosomes with your own server system.

After imputation, the following command can be used to concatenate all single-cell imputed matrices. You need to provide a list of imputed files need to be concatenated.
```
ls {impute_dir}/*{imputation_mode}_{chromosome}.hdf5 > {impute_file_list}
hicluster embed-concatcell-chr --cell_list {impute_file_list} --outprefix {embed_dir}{imputation_mode}_{chromosome} --res ${resolution}

ls {embed_dir}{imputation_mode}_*npy > {embed_file_list}
hicluster embed-mergechr --embed_list {embed_file_list} --outprefix {embed_dir}{imputation_mode}
```
where the imputation_mode is a name that can be defined by the user in the impute-cell command, and default to be 'pad?_std?_rp?_sqrtvc' based on the imputation parameters.

The embedding generated here could be further used for batch-effect correction, clustering, and visulization.

### Loop calling

The loop calling framework is modified from SnapHiC. When using these functions, please cite [this work](https://www.biorxiv.org/content/10.1101/2020.12.13.422543v1). The algorithm uses single-cell imputed matrices as input. It first normalize the matrices by distance decay, and compute the difference between each normalized element and its local background next. Finally the single-cell level matrices of the same cell type are merged, and paired t-test across cells and other ad-hoc filters are applied to identify loops.  



### Merge single cells

When you have a list of cells that you're interested in, either determined from clustering result or from expriment labels, you can merge the contact matrices of all those single cells to generate a pseudo bulk data.
If you want to merge the scHiCluster imputed matrices, the command should be
```
Q = merge_gpu(network, c, res, pad=1, rp=0.5, prct=-1)
```
If you want to merge the raw contact matrices directly, you can use
```
Q = merge_raw(network, c, res)
```
c represents the chromosome, and res represent the resolution.

### Single cell domain calling

We haven't incooperate the domain calling function directly in the package, but you can output the imputed contact matrices or raw contact matrices and use them as input for other domain calling software. Here we provided two output format, including a sparse matrix format, and the format that Topdom required. You can use
```
output_topdom(cell, c, Q, res)
output_sparse(cell, c, Q, res)
```
where cell is the path and cell name of the output file. And the function can be used for both single cell matrices and merged matrices. One may also need to notice that 
1. We suggest to use only cells with greater than 50k contacts for this analysis.
2. The imputed matrices are usually not sparse even for a single cell, so you may consider save with less float points, or only save the contacts within certain distance from the diagnol when the matrices are large.

After we have a contact matrix in text format, we can run Topdom(http://zhoulab.usc.edu/TopDom/) to find domains. Run
```
source('TopDom_v0.0.2.R')
tad = TopDom(matrix.file = 'cell_1_chr1.matrix', window.size = 5)
write.table(tad$bed[1:dim(tad$bed)[1],2:4], file='cell_1_chr1.w5.domain', quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
```
in R to get the domain calling result. The output format will be 

> 0 325000  gap
> 325000 425000 domain
> 425000 650000 domain

### Differential domain boundaries between cell types

Suppose you have the domain calling results output by Topdom in each single cells, you can identify the differential domain boundaries across cell types on a single chromosome by  
```
sc_dom, dom_prob, bins, pvalue = diff_dom(domlist, cluster, celltypelist, res, chrom, chromsize)
```
Input:
domlist: a list of all the topdom output.
cluster: the cluster assignment of each cell with the same order of domlist.
celltypelist: a list of all possible cell type label.
res: the resolution of domains.
chrom: the chromosome.

Returns:
sc_dom: domain boundaries in each single cell indicated by a binary matrix.
dom_prob: domain boundary frequency in each cluster with the same order as in celltypelist.
bins: The tested bins, since the bins with 0.0 or 1.0 domain frequency in any of the clusters were not tested.
pvalue: The p-value of each bin in bins.

We define the function on single chromosome level to allow user to parallel processing all the chromosomes together. After getting results from all chromosomes, FDR correction and other thresholding on prob_dom should be applied. 
```
from statsmodels.sandbox.stats.multicomp import multipletests as FDR
fdr = FDR(pvalue, 0.01, 'fdr_bh')[1]
sigfilter = filter_bins(dom_prob, fdr, fdr_cutoff, diff_cutoff, max_cutoff, min_cutoff)
```
diff_cutoff is the minimal boundary probability difference between the cluster with highest boundary probability minus the cluster with lowest boundary probability.
max_cutoff is the minimal boundary probability of the cluster with highest boundary probability.
min_cutoff is the maximal boundary probability of the cluster with lowest boundary probability.

Since the domain boundaries sometimes shift for a few bins, we suggest to filter the differential domains considering the flanking bins. You can use the following function to compute the p-value again using the flanking w bins including the bin being tested.
```
flank_dom_prob, flank_pvalue, bins = diff_dom_flank(sigbins, scdom, cluster, celltypelist, res, chrom, chromsize, w)
```
After concatenating result from all chromosomes. You can run
```
sigfilter = filter_bins(flank_dom_prob, fdr, fdr_cutoff, diff_cutoff, max_cutoff, min_cutoff)
sigbins = bins[sigfilter]
```
to get the final differential domain boundary calling.
