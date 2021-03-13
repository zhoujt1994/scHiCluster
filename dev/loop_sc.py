import sys
import h5py
import time
import cv2
cv2.useOptimized()
import argparse
import numpy as np
from scipy.sparse import load_npz, save_npz
from sklearn.preprocessing import RobustScaler

parser = argparse.ArgumentParser()
parser.add_argument('--outdir', type=str, default=None, help='Directory of imputed matrix')
parser.add_argument('--cell', type=str, default=None, help='Specific identifier of a cell')
parser.add_argument('--chrom', type=str, default=None, help='Chromosome imputed')
parser.add_argument('--impute_mode', type=str, default=None, help='Suffix of imputed matrix file names')
parser.add_argument('--res', type=int, default=None, help='Bin size as integer to generate contact matrix')
parser.add_argument('--dist', type=int, default=10000000, help='Maximum distance threshold of contacts to use')
parser.add_argument('--cap', type=int, default=5, help='Trim Z-scores over the threshold')
parser.add_argument('--pad', type=int, default=5, help='One direction size of larger square for donut background')
parser.add_argument('--gap', type=int, default=2, help='One direction size of smaller square for donut background')
parser.add_argument('--norm_mode', type=str, default='dist_trim', help='Suffix of normalized file names')
opt = parser.parse_args()

if opt.chrom[:3]=='chr':
    c = opt.chrom[3:]
else:
    c = opt.chrom

E = load_npz(opt.outdir + opt.cell + '_chr' + c + '_' + opt.impute_mode + '.npz')
start_time = time.time()
ave, std, top, count = np.zeros((4, opt.dist // opt.res + opt.pad + 1))
for i in range(opt.dist // opt.res + opt.pad + 1):
	tmp = E.diagonal(i)
	top[i] = np.percentile(tmp,99)
	tmp[tmp>top[i]] = top[i]
	ave[i] = np.mean(tmp)
	std[i] = np.std(tmp)
	count[i] = np.sum(tmp>0)

print('Curve', time.time() - start_time, '#Nonzero', np.sum(count))

start_time = time.time()
E.data = np.min([E.data, top[E.col - E.row]], axis=0)
E = E.astype(np.float32).toarray()
tmp = E[idx]
tmp = (tmp - ave[idx[1] - idx[0]]) / std[idx[1] - idx[0]]
tmp[count[idx[1] - idx[0]] < 100] = 0
tmp[std[idx[1] - idx[0]]==0] = 0
tmp[tmp > opt.cap] = opt.cap
tmp[tmp < -opt.cap] = -opt.cap
E[idx] = tmp.copy()
print('Norm', time.time() - start_time)

start_time = time.time()
w = opt.pad * 2 + 1
kernel = np.ones((w, w), np.float32)
kernel[(opt.pad - opt.gap):(opt.pad + opt.gap + 1), (opt.pad - opt.gap):(opt.pad + opt.gap + 1)] = 0
kernel = kernel / np.sum(kernel)
N = cv2.filter2D(E, -1, kernel=kernel)
E = csr_matrix(np.around(E, decimals=6))
N = csr_matrix(np.around(N, decimals=6)).multiply(mask)
N = E - N
print('Bkg', time.time() - start_time)

save_npz(opt.outdir + opt.cell + '_chr' + c + '_' + opt.impute_mode + '_' + opt.norm_mode + '.E.npz', E)
save_npz(opt.outdir + opt.cell + '_chr' + c + '_' + opt.impute_mode + '_' + opt.norm_mode + '.T.npz', N)


