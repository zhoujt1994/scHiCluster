# command time ~/miniconda3/envs/r_env/bin/Rscript /gale/ddn/snm3C/humanPFC/code/domain_topdom_cell.R ${sample} ${c} pad2_std1_rp0.5_sqrtvc 10 /gale/ddn/snm3C/humanPFC/smoothed_matrix/${res0}b_resolution/ /gale/ddn/snm3C/humanPFC/smoothed_matrix/${res0}b_resolution/

library(rhdf5)
library(Matrix)
library(data.table)
source('/gale/netapp/home/zhoujt/software/TopDom_v0.0.2.R')

args = commandArgs(trailingOnly=TRUE)
cell = args[1]
if (args[2][1:3]=='chr'){chrom = args[2]}else{chrom = paste('chr', args[2], sep='')}
mode = args[3]
ws = args[4]
indir = args[5]
outdir = args[6]
fmat = paste(indir, chrom, '/', cell, '_', chrom, '_', mode, '.hdf5', sep = '')
fbin = paste(indir, 'bins/', chrom, '.bed', sep='')
bins = fread(fbin)
bins = data.frame(bins)
colnames(bins) <- c("chr", "from.coord", "to.coord")
n_bins = nrow(bins)
indices = as.numeric(h5read(fmat, 'Matrix/indices'))
indptr = as.numeric(h5read(fmat, 'Matrix/indptr'))
count = as.numeric(h5read(fmat, 'Matrix/data'))
matrix.data = as.matrix(sparseMatrix(j=indices+1, p=indptr, x=count, dims=c(n_bins, n_bins)))
tad = TopDom(matrix.data, bins, window.size=strtoi(ws))
fout = paste(outdir, chrom, '/', cell, '_', chrom, '_', mode, '.w', ws, '.domain.bed', sep = '')
write.table(tad$bed, file=fout, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
