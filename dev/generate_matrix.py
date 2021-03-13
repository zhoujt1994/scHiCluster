import argparse
import numpy as np
from collections import Counter

parser = argparse.ArgumentParser()
parser.add_argument('--infile', type=str, default=None, help='Path to the short format contact file')
parser.add_argument('--outdir', type=str, default=None, help='Output directory end with /')
parser.add_argument('--cell', type=str, default=None, help='Specific identifier of a cell')
parser.add_argument('--res', type=int, default=None, help='Bin size as integer to generate contact matrix')
parser.add_argument('--genome', type=str, default=None, help='Genome assembly version') ## mm10 or hg38
parser.add_argument('--dist', type=int, default=None, help='Minimum distance threshold of contacts to use')
opt = parser.parse_args()

fin = open(opt.infile, 'r')
if opt.res>=1000000:
	k = str(opt.res/1000000)+'mb'
elif res>=1000:
	k = str(opt.res/1000)+'kb'

if opt.genome[:2]=='hg':
	chrom = [str(i+1) for i in range(22)]
elif opt.genome[:2]=='mm':
	chrom = [str(i+1) for i in range(19)]

count = Counter()
for line in fin:
	tmp = line.strip().split('\t')
	tmp[1], tmp[5] = tmp[1][3:], tmp[5][3:]
	if (np.abs(int(tmp[6])-int(tmp[2])) >= opt.dist) and (tmp[1]==tmp[5]) and (tmp[1] in chrom):
		tmp[2], tmp[6] = int(tmp[2]) // res, int(tmp[6]) // res
		if tmp[2]>tmp[6]:
			tmp[2], tmp[6] = tmp[6], tmp[2]
		if tmp[2]!=tmp[6]:
			pos = '-'.join([tmp[1], str(tmp[2]), str(tmp[6])])
			count[pos] += 1

fin.close()
fout = {c:open(opt.outdir + 'chr' + c + '/' + opt.cell + '_chr' + c + '.txt', 'w') for c in chrom}
for key in count:
	tmp = key.split('-')
	fout[tmp[0]].write('{0}\t{1}\t{2}.0\n'.format(tmp[1], tmp[2], count[key]))

for c in chrom:
	fout[c].close()


