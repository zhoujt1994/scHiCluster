library(rhdf5)
library(Matrix)
library(data.table)
source('TopDom_v0.0.2.R')

args = commandArgs(trailingOnly=TRUE)
cell = args[1]
chrom = paste('chr', args[2], sep='')
ws = args[3]
indir = ''
outdir = ''
fmat = paste(indir, chrom, '/', cell, '_', chrom, '.hdf5', sep = '')
fbin = paste(indir, 'bins/', chrom, '.bed', sep='')
bins = fread(fbin)
bins = data.frame(bins)
colnames(bins) <- c("chr", "from.coord", "to.coord")
n_bins = nrow(bins)
idx = as.numeric(h5read(fmat, 'row'))
idy = as.numeric(h5read(fmat, 'col'))
count = as.numeric(h5read(fmat, 'data'))
matrix.data = as.matrix(sparseMatrix(i=idx+1, j=idy+1, x=count, dims=c(n_bins, n_bins)))
tad = TopDom(matrix.data, bins, window.size=strtoi(ws))
fout = paste(outdir, chrom, '/', cell, '_', chrom, '.w', ws, '.domain', sep = '')
write.table(tad$bed, file=fout, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
