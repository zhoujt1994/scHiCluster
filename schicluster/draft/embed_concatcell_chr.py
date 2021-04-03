# command time python /gale/ddn/snm3C/humanPFC/code/embed_concatcell_chr.py --cell_list /gale/ddn/snm3C/humanPFC/smoothed_matrix/filelist/imputelist_pad1_std1_rp0.5_sqrtvc_chr${c}.txt --outprefix /gale/ddn/snm3C/humanPFC/smoothed_matrix/${res0}b_resolution/merged/pad1_std1_rp0.5_sqrtvc_chr${c} --chrom ${c} --res ${res}

import time
import h5py
import argparse
import numpy as np
from scipy.sparse import load_npz, save_npz, csr_matrix, vstack
from sklearn.decomposition import TruncatedSVD

def embed_concatcell_chr(cell_list, outprefix, res, dist=10000000, save_raw=True, dim=50):

	celllist = np.loadtxt(cell_list, dtype=np.str)

	with h5py.File(celllist[0], 'r') as f:
		ngene = f['Matrix'].attrs['shape'][0]
	idx = np.triu_indices(ngene, k=1)
	idxfilter = np.array([(yy - xx) < (dist / res + 1) for xx,yy in zip(idx[0], idx[1])])
	idx = (idx[0][idxfilter], idx[1][idxfilter])

	start_time = time.time()
	# matrix = np.zeros((len(celllist), np.sum(idxfilter)))
	matrix = []
	for i,cell in enumerate(celllist):
		with h5py.File(cell, 'r') as f:
			g = f['Matrix']
			A = csr_matrix((g['data'][()], g['indices'][()], g['indptr'][()]), g.attrs['shape'])
		# matrix[i] = A[idx]
		matrix.append(csr_matrix(A[idx]))
		if i%100==0:
			print(i, 'cells loaded', time.time() - start_time, 'seconds')

	matrix = vstack(matrix)

	if save_raw:
		save_npz(f'{outprefix}.npz', matrix)

	scalefactor = 100000
	matrix.data = matrix.data * scalefactor
	svd = TruncatedSVD(n_components=dim, algorithm='arpack')
	matrix_reduce = svd.fit_transform(matrix)
	matrix_reduce = matrix_reduce / svd.singular_values_
	np.save(f'{outprefix}.svd{dim}.npy', matrix_reduce)
	return
'''
parser = argparse.ArgumentParser()
parser.add_argument('--cell_list', type=str, default=None, help='Full path of a file containing the full path to all imputed matrices to be concatenate')
parser.add_argument('--outprefix', type=str, default=None, help='Prefix of concatenated matrix including directory')
parser.add_argument('--res', type=int, default=None, help='Bin size as integer to generate contact matrix')
parser.add_argument('--dist', type=int, default=10000000, help='Maximum distance threshold of contacts to use')
parser.add_argument('--skip_raw', dest='save_raw', action='store_false', help='Not to save cell-by-feature matrix before SVD')
parser.set_defaults(save_raw=True)
parser.add_argument('--dim', type=int, default=50, help='Number of dimensions to return from SVD')
opt = parser.parse_args()

concat_cell(opt.cell_list, opt.outprefix, opt.res, opt.dist, opt.save_raw, opt.dim)
'''
