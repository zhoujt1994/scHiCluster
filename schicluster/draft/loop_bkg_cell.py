# command time python /gale/ddn/snm3C/humanPFC/code/loop_bkg_cell.py --indir /gale/ddn/snm3C/humanPFC/smoothed_matrix/${res0}b_resolution/ --cell ${sample} --chrom ${c} --res ${res} --impute_mode pad1_std1_rp0.5_sqrtvc

import sys
import h5py
import time
import cv2
cv2.useOptimized()
import argparse
import numpy as np
from scipy.sparse import save_npz, csr_matrix
from sklearn.preprocessing import RobustScaler

def loop_bkg_cell(indir, cell, chrom, impute_mode, res, 
			dist=10050000, cap=5, pad=5, gap=2, norm_mode='dist_trim'):

	if chrom[:3]=='chr':
		c = chrom[3:]
	else:
		c = chrom

	start_time = time.time()
	with h5py.File(f'{indir}chr{c}/{cell}_chr{c}_{impute_mode}.hdf5', 'r') as f:
		g = f['Matrix']
		E = csr_matrix((g['data'][()], g['indices'][()], g['indptr'][()]), g.attrs['shape']).tocoo()

	print('Load', time.time() - start_time)

	# TODO numba optimize
	start_time = time.time()
	ave, std, top, count = np.zeros((4, dist // res + pad + 1))
	for i in range(dist // res + pad + 1):
		tmp = E.diagonal(i)
		top[i] = np.percentile(tmp,99)
		tmp[tmp>top[i]] = top[i]
		ave[i] = np.mean(tmp)
		std[i] = np.std(tmp)
		count[i] = np.sum(tmp>0)

	print('Curve', time.time() - start_time, '#Nonzero', np.sum(count))

	idx = np.triu_indices(E.shape[0], 0)
	idxfilter = ((idx[1] - idx[0]) < (dist // res + 1))
	idx = (idx[0][idxfilter], idx[1][idxfilter])
	mask = csr_matrix((np.ones(len(idx[0])), (idx[0], idx[1])), E.shape)

	start_time = time.time()
	E.data = np.min([E.data, top[E.col - E.row]], axis=0)
	E = E.astype(np.float32).toarray()
	tmp = E[idx]
	tmp = (tmp - ave[idx[1] - idx[0]]) / std[idx[1] - idx[0]]
	tmp[count[idx[1] - idx[0]] < 100] = 0
	tmp[std[idx[1] - idx[0]]==0] = 0
	tmp[tmp > cap] = cap
	tmp[tmp < -cap] = -cap
	E[idx] = tmp.copy()
	print('Norm', time.time() - start_time)

	start_time = time.time()
	w = pad * 2 + 1
	kernel = np.ones((w, w), np.float32)
	kernel[(pad - gap):(pad + gap + 1), (pad - gap):(pad + gap + 1)] = 0
	kernel = kernel / np.sum(kernel)
	N = cv2.filter2D(E, -1, kernel=kernel)
	E = csr_matrix(np.around(E, decimals=6))
	N = csr_matrix(np.around(N, decimals=6)).multiply(mask)
	N = E - N
	print('Bkg', time.time() - start_time)

	save_npz(f'{indir}chr{c}/{cell}_chr{c}_{impute_mode}_{norm_mode}.E.npz', E)
	save_npz(f'{indir}chr{c}/{cell}_chr{c}_{impute_mode}_{norm_mode}.T.npz', N)
	return

parser = argparse.ArgumentParser()
parser.add_argument('--indir', type=str, default=None, help='Directory of imputed matrix')
parser.add_argument('--cell', type=str, default=None, help='Specific identifier of a cell')
parser.add_argument('--chrom', type=str, default=None, help='Chromosome imputed')
parser.add_argument('--impute_mode', type=str, default=None, help='Suffix of imputed matrix file names')
parser.add_argument('--res', type=int, default=None, help='Bin size as integer to generate contact matrix')

parser.add_argument('--dist', type=int, default=10050000, help='Maximum distance threshold of contacts to use')
parser.add_argument('--cap', type=int, default=5, help='Trim Z-scores over the threshold')
parser.add_argument('--pad', type=int, default=5, help='One direction size of larger square for donut background')
parser.add_argument('--gap', type=int, default=2, help='One direction size of smaller square for donut background')
parser.add_argument('--norm_mode', type=str, default='dist_trim', help='Suffix of normalized file names')
opt = parser.parse_args()

loop_bkg_cell(opt.indir, opt.cell, opt.chrom, opt.impute_mode, opt.res, 
		opt.dist, opt.cap, opt.pad, opt.gap, opt.norm_mode)
