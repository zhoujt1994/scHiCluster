# python /gale/ddn/snm3C/humanPFC/code/generatematrix_cell.py --infile /gale/netapp/seq4/illumina_runs/ENTEx/Pool_BX_BY/mapping/output/TrCo_A9HOW_Plate1-1-H22/hic/TrCo_A9HOW_Plate1-1-H22-A13.3C.sorted_contacts_old.txt --outdir /gale/ddn/snm3C/humanPFC/cell_matrix/test/10kb_resolution/ --chrom_file /gale/netapp/home/zhoujt/genome/hg19/hg19.autosomal.chrom.sizes --split_file /gale/netapp/home/zhoujt/genome/hg19/hg19.chrsplit.bed --res 10000 --cell TrCo_A9HOW_Plate1-1-H22-A13
# python /gale/ddn/snm3C/humanPFC/code/generatematrix_cell.py --infile /gale/netapp/seq4/illumina_runs/ENTEx/Pool_BX_BY/mapping/output/TrCo_A9HOW_Plate1-1-H22/hic/TrCo_A9HOW_Plate1-1-H22-A13.3C.sorted_contacts_old.txt --outdir /gale/ddn/snm3C/humanPFC/cell_matrix/test/100kb_resolution/ --chrom_file /gale/netapp/home/zhoujt/genome/hg19/hg19.autosomal.chrom.sizes --res 100000 --cell TrCo_A9HOW_Plate1-1-H22-A13

import os
import argparse
import numpy as np
import pandas as pd
#from collections import Counter

def generatematrix_cell(infile, outdir, cell, res, chrom_file, 
					chr1=1, pos1=2, chr2=5, pos2=6, split_file=None, dist=2500):

	chrom = np.loadtxt(chrom_file, dtype=np.str)[:,0]

	# add p/q arm to split chromosomes
	if not split_file:
		chrom_split = chrom.copy() 
	else:
		splitbed = pd.read_csv(split_file, sep='\t', header=None, index_col=0)
		chrom_split = np.concatenate([[c+'p', c+'q'] if c in splitbed.index else [c] for c in chrom])

	# create a folder per chromosome, containing all cells
	for c in chrom_split:
		if c[:3]=='chr':
			os.makedirs(f'{outdir}{c}/', exist_ok=True)
		else:
			os.makedirs(f'{outdir}chr{c}/', exist_ok=True)

	data = pd.read_csv(infile, sep='\t', header=None)
	# make pos1 < pos2
	data.loc[data[pos1]>data[pos2], [pos1, pos2]] = data.loc[data[pos1]>data[pos2], [pos2, pos1]]
	# filter intra autosomal contacts with insertion size threshold
	data = data[(data[chr1]==data[chr2]) & (data[chr1].isin(chrom)) & (data[pos2] - data[pos1] > dist)]
	if split_file:
		splitfilter = data[chr1].isin(splitbed.index)
		datasplit = data[splitfilter]
		# end position of p and start position of q
		datasplit['p'] = splitbed.loc[datasplit[chr1], 1].values
		datasplit['q'] = splitbed.loc[datasplit[chr1], 2].values
		datasplit.loc[datasplit[pos1]<datasplit['p'], chr1] += 'p'
		datasplit.loc[datasplit[pos2]<datasplit['p'], chr2] += 'p'	
		datasplit.loc[datasplit[pos1]>datasplit['q'], chr1] += 'q'
		datasplit.loc[datasplit[pos2]>datasplit['q'], chr2] += 'q'
		# fully divide to ensure bins consistency between p and q
		datasplit.loc[datasplit[pos1]>datasplit['q'], pos1] -= datasplit['q'] // res * res
		datasplit.loc[datasplit[pos2]>datasplit['q'], pos2] -= datasplit['q'] // res * res
		data = pd.concat([datasplit[data.columns], data[~splitfilter]], axis=0)
	data[[pos1, pos2]] = data[[pos1,pos2]] // res
	data = data[(data[chr1]==data[chr2]) & (data[chr1].isin(chrom_split)) & (data[pos1]!=data[pos2])]
	# count contacts at each pixel
	data = data.groupby(by=[chr1, pos1, pos2])[chr2].count().reset_index()
	for c, tmp in data.groupby(by=chr1):
		if c[:3]=='chr':
			tmp[[pos1, pos2, chr2]].astype(int).to_csv(f'{outdir}{c}/{cell}_{c}.txt', sep='\t', header=False, index=False)
		else:
			tmp[[pos1, pos2, chr2]].astype(int).to_csv(f'{outdir}{c}/{cell}_chr{c}.txt', sep='\t', header=False, index=False)
	'''
	count = Counter()
	with open(infile, 'r') as fin:
		for line in fin:
			tmp = line.strip().split('\t')
			c1, c2, p1, p2 = tmp[chr1], tmp[chr2], int(tmp[pos1]), int(tmp[pos2])
			if p1>p2:
				p1, p2 = p2, p1
			if split_file:
				if c1 in splitbed.index:
					if p1 < splitbed.loc[c1,0]:
						c1 = c1 + 'p'
					elif p1 > splitbed.loc[c1,1]:
						c1 = c1 + 'q'
						p1 = p1 - splitbed.loc[c1,1]
				if c2 in splitbed.index:
					if p2 < splitbed.loc[c2,0]:
						c2 = c2 + 'p'
					elif p2 > splitbed.loc[c2,1]:
						c2 = c2 + 'q'
						p2 = p2 - splitbed.loc[c2,1]
			if (c1==c2) and (c1 in chrom_split) and (p2 - p1 >= dist):
				p1, p2 = p1 // res, p2 // res
				if p1!=p2:
					pos = '-'.join([c1, str(p1), str(p2)])
					count[pos] += 1
	fout = {c:open(outdir + 'chr' + c + '/' + cell + '_chr' + c + '.txt', 'w') for c in chrom_split}
	for key in count:
		tmp = key.split('-')
		fout[tmp[0]].write('{0}\t{1}\t{2}.0\n'.format(tmp[1], tmp[2], count[key]))
	for c in chrom:
		fout[c].close()
	'''
	return

parser = argparse.ArgumentParser()
parser.add_argument('--infile', type=str, default=None, help='Path to the short format contact file')
parser.add_argument('--chr1', type=int, default=2, help='Column of fragment1 chromosome')
parser.add_argument('--pos1', type=int, default=3, help='Column of fragment1 position')
parser.add_argument('--chr2', type=int, default=6, help='Column of fragment2 chromosome')
parser.add_argument('--pos2', type=int, default=7, help='Column of fragment2 position')
parser.add_argument('--outdir', type=str, default=None, help='Output directory end with /')
parser.add_argument('--cell', type=str, default=None, help='Specific identifier of a cell')
parser.add_argument('--res', type=int, default=None, help='Bin size as integer to generate contact matrix')
parser.add_argument('--chrom_file', type=str, default=None, help='Path to the chromosome size files containing all chromosome to be analyzed as the first column') ## mm10 or hg38
parser.add_argument('--split_file', type=str, default=None, help='Path to the bed file containing all chromosomes need to split and one region per chromosome as splitting point')
parser.add_argument('--dist', type=int, default=2500, help='Minimum distance threshold of contacts to use')
opt = parser.parse_args()

generatematrix_cell(opt.infile, opt.outdir, opt.cell, opt.res, opt.chrom_file, 
				opt.chr1-1, opt.pos1-1, opt.chr2-1, opt.pos2-1, opt.split_file, opt.dist)
