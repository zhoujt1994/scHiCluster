# command time python /gale/ddn/snm3C/humanPFC/code/concat_cell.py --indir /gale/ddn/snm3C/humanPFC/smoothed_matrix/${res0}b_resolution/ --cell_list /gale/ddn/snm3C/humanPFC/smoothed_matrix/celllist_long.txt --chrom ${SGE_TASK_ID} --mode pad1_std1_rp0.5_sqrtvc --res ${res}

import time
import h5py
import argparse
import numpy as np
from scipy.sparse import load_npz, save_npz, csr_matrix, vstack
from sklearn.decomposition import TruncatedSVD

def concat_cell(indir, cell_list, chrom, res, mode, dist=10000000, save_raw=True):
	if chrom[:3]=='chr':
	    c = chrom[3:]
	else:
	    c = chrom

	celllist = np.loadtxt(cell_list, dtype=np.str)

	with h5py.File(f'{indir}chr{c}/{celllist[0]}_chr{c}_{mode}.hdf5', 'r') as f:
		ngene = f['Matrix'].attrs['shape'][0]
	idx = np.triu_indices(ngene, k=1)
	idxfilter = np.array([(yy - xx) < (dist / res + 1) for xx,yy in zip(idx[0], idx[1])])
	idx = (idx[0][idxfilter], idx[1][idxfilter])

	start_time = time.time()
	matrix = np.zeros((len(celllist), np.sum(idxfilter)))
	for i,cell in enumerate(celllist):
		with h5py.File(f'{indir}chr{c}/{cell}_chr{c}_{mode}.hdf5', 'r') as f:
			g = f['Matrix']
			A = csr_matrix((g['data'][()], g['indices'][()], g['indptr'][()]), g.attrs['shape'])
		matrix[i] = csr_matrix(A[idx])
		if i%100==0:
			print(i, 'cells loaded', time.time() - start_time, 'seconds')

	matrix = vstack(matrix)

	if save_raw:
		save_npz(f'{indir}merged/{mode}_chr{c}.npz', matrix)

	scalefactor = 100000
	matrix.data = matrix.data * scalefactor
	svd = TruncatedSVD(n_components=50, algorithm='arpack')
	matrix_reduce = svd.fit_transform(matrix)
	matrix_reduce = matrix_reduce / svd.singular_values_
	np.save(f'{indir}merged/{mode}_chr{c}.svd50.npy', matrix_reduce)

parser = argparse.ArgumentParser()
parser.add_argument('--indir', type=str, default=None, help='Directory of imputed matrices end with /')
parser.add_argument('--cell_list', type=str, default=None, help='Full path of a file containing a list of cell identifiers to be concatenate')
parser.add_argument('--chrom', type=str, default=None, help='Chromosome to impute')
parser.add_argument('--mode', type=str, default=None, help='Suffix of imputed matrix file names')
parser.add_argument('--res', type=int, default=None, help='Bin size as integer to generate contact matrix')
parser.add_argument('--dist', type=int, default=10000000, help='Maximum distance threshold of contacts to use')
parser.add_argument('--save_raw', type=bool, default=True, help='Whether to save cell-by-feature matrix before SVD')
opt = parser.parse_args()

concat_cell(opt.indir, opt.cell_list, opt.chrom, opt.res, opt.mode, opt.dist, opt.save_raw)
