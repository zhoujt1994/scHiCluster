import sys
import pandas as pd

chrom_file = sys.argv[4]
chromsize = pd.read_csv(chrom_file, sep='\t', header=None, index_col=0).to_dict()[1]
res = 500000
size = [0]
for c in chromsize:
	size.append(int(chromsize[c] // res) + 1 + size[-1])

size = {xx:yy for xx,yy in zip(chromsize, size)}

infile = sys.argv[1]
outdir = sys.argv[2]
cell = sys.argv[3]
data = pd.read_csv(infile, sep='\t', header=None, index_col=None)
data = data[data[4]==data[5]]
data[4] = data[4].str.split('_', expand=True)[1]
data = data[data[4].isin(chromsize.keys())]
data[5] = data[4].map(size)
data[[0,1]] = data[[0,1]].values - data[5].values[:,None]
for c, tmp in data.groupby(by=4):
	tmp[[0,1,2]].astype(int).to_csv(f'{outdir}{c}/{cell}_{c}.txt', sep='\t', header=False, index=False)

tmp = data[4].unique()
for c in chromsize:
	if not (c in tmp):
		open(f'{outdir}{c}/{cell}_{c}.txt', 'w').close()
