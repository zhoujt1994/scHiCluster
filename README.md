# scHiCluster
## Installation
Running scHiCluster requries numpy, scipy, scikit-learn.
Running the GPU version requires pytorch in addition.
First, creat a new conda environment and activate it by
```
conda create --name schicluster python==3.6.1
source activate schicluster
```
Then install the prerequsite packages by
```
conda install -c anaconda numpy scipy scikit-learn
```
And install pytorch for your cuda and python version according to the instruction on https://pytorch.org/get-started/locally/
Finally install scHiCluster by
```
git clone https://github.com/zhoujt1994/HiCluster.git
cd HiCluster
pip install .
```

## Usage
### Clustering
#### Data format
You need to prepare your Hi-C contact matrices in a sparse matrix format. You need three columns separated by tab, representing the interacting bins and the number of reads supporting the interaction. The name of the file need to be your cell name, followed by the chromosome name.
For example, at 1mb resolution, if there are 10 reads supporting the interaction between chr1:1000000-2000000 and chr1:5000000-6000000, then you should have a line as
> 1 5 10.0

in a file named as cellname_chr1.txt

Then you will need a sample list providing to the program that contains all the cell names you want to analyze, without the chromosome names. For instance, you need to provide an array variable named network, where each element is format like
> directory/cell_i

the program will load the files directory/cell_i_chr1.txt to cell_i_chr22.txt for human cells.
And you also need a dictionary variable chromsize, to provide the length of each chromosome in the samples you want to analyze.

With those information ready, scHiCluster can be run by
```
cluster, embedding = hicluster_gpu(network, chromsize, nc=nc, pad=1, rp=0.5, prct=20)
```
or
```
cluster, embedding = hicluster_cpu(network, chromsize, nc=nc, pad=1, rp=0.5, prct=20, ncpus=5)
```
nc represents the number of clusters
pad represent the window size for linear convolution
rp represent the 1 - restart probablity
prct represent the percentage of largest values that kept for each cell

The function will return the discrete cluster assignment and the embedding, which can be used for other analysis and visualization.

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
