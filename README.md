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

When having the merged matrix, you can output it in a sparse matrix format, or the format that Topdom required as input by
```
output_topdom(cell, c, Q, res)
output_sparse(cell, c, Q, res)
```
where cell is the path and cell name of the output file.
