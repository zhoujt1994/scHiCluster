# command time python /gale/ddn/snm3C/humanPFC/code/domain_insulation_cell.py --indir /gale/ddn/snm3C/humanPFC/smoothed_matrix/${res0}b_resolution/ --cell ${sample} --chrom ${c} --mode pad2_std1_rp0.5_sqrtvc --w 10

import h5py
import argparse
import numpy as np
from scipy.sparse import csr_matrix

def domain_insulation_cell(indir, cell, chrom, mode, w=10):

	if chrom[:3]=='chr':
		c = chrom
	else:
		c = 'chr' + chrom
	with h5py.File(f'{indir}{c}/{cell}_{c}_{mode}.hdf5', 'r') as f:
		g = f['Matrix']
		A = csr_matrix((g['data'][()], g['indices'][()], g['indptr'][()]), g.attrs['shape'])
	score = np.ones(A.shape[0])
	for i in range(1,A.shape[0]):
		if i<w:
			intra = (A[:i, :i].sum() + A[i:(i+w), i:(i+w)].sum()) / (i*(i+1)/2 + (w*(w+1)/2))
			inter = A[:i, i:(i+w)].sum() / (i * (i+w))
		else:
			intra = (A[(i-w):i, (i-w):i].sum() + A[i:(i+w), i:(i+w)].sum()) / (w*(w+1))
			inter = A[(i-w):i, i:(i+w)].sum() / (w*w)
		score[i] = inter / (inter + intra)
	np.save(f'{indir}{c}/{cell}_{c}_{mode}.w{w}.ins.npy', score)
	return

parser = argparse.ArgumentParser()
parser.add_argument('--indir', type=str, default=None, help='Directory of imputed matrices end with /')
parser.add_argument('--cell', type=str, default=None, help='Full path of a file containing a list of cell identifiers to be concatenate')
parser.add_argument('--chrom', type=str, default=None, help='Chromosome to impute')
parser.add_argument('--mode', type=str, default=None, help='Suffix of imputed matrix file names')
parser.add_argument('--w', type=int, default=10, help='Window size for insulation score')
opt = parser.parse_args()

domain_insulation_cell(opt.indir, opt.cell, opt.chrom, opt.mode, opt.w)
