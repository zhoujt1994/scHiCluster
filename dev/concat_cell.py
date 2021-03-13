import time
import h5py
import argparse
import numpy as np
from scipy.sparse import load_npz, save_npz, csr_matrix, vstack
from sklearn.decomposition import TruncatedSVD

parser = argparse.ArgumentParser()
parser.add_argument('--indir', type=str, default=None, help='Directory of imputed matrices end with /')
parser.add_argument('--chrom', type=str, default=None, help='Chromosome to impute')
parser.add_argument('--mode', type=str, default=None, help='Suffix of imputed matrix file names')
parser.add_argument('--res', type=int, default=None, help='Bin size as integer to generate contact matrix')
parser.add_argument('--dist', type=int, default=10000000, help='Maximum distance threshold of contacts to use')
parser.add_argument('--save_raw', type=bool, default=True, help='Whether to save cell-by-feature matrix before SVD')
opt = parser.parse_args()

if opt.chrom[:3]=='chr':
    c = opt.chrom[3:]
else:
    c = opt.chrom

celllist = np.loadtxt(opt.indir + '../cell_list.txt', dtype=np.str)

f = h5py.File(opt.indir + 'chr' + c + '/' + celllist[0] + '_chr' + c + '_' + opt.mode + '.hdf5', 'r')
ngene = f['Matrix'].attrs['shape'][0]
f.close()
idx = np.triu_indices(ngene, k=1)
idxfilter = np.array([(yy - xx) < (opt.dist / opt.res + 1) for xx,yy in zip(idx[0], idx[1])])
idx = (idx[0][idxfilter], idx[1][idxfilter])

start_time = time.time()
matrix = []
for i,cell in enumerate(celllist):
	f = h5py.File(opt.indir + 'chr' + c + '/' + cell + '_chr' + c + '_' + opt.mode + '.hdf5', 'r')
	g = f['Matrix']
	A = csr_matrix((g['data'][()], g['indices'][()], f['indptr'][()]), g.attrs['shape'])
	f.close()
	matrix.append(csr_matrix(A[idx]))
	if i%100==0:
		print(i, 'cells loaded', end_time-start_time, 'seconds')

matrix = vstack(matrix)

if opt.save_raw:
	save_npz(opt.indir + 'merged/' + opt.mode + '_chr' + c + '.npz', matrix)

scalefactor = 100000
matrix.data = matrix.data * scalefactor
svd = TruncatedSVD(n_components=50, algorithm='arpack')
matrix_reduce = svd.fit_transform(matrix)
matrix_reduce = matrix_reduce / svd.singular_values_
np.save(opt.indir + 'merged/' + opt.mode + '_chr' + c + '.svd50.npy', matrix_reduce)
