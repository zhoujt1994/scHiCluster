import time
import cv2
import argparse
import numpy as np
import pandas as pd
from scipy import stats
from scipy.sparse import save_npz, load_npz, csr_matrix, coo_matrix

parser = argparse.ArgumentParser()
parser.add_argument('--indir', type=str, default=None, help='Directory of imputed matrices end with /')
parser.add_argument('--cell_list', type=str, default=None, help='List of cell identifiers to be merged')
parser.add_argument('--group', type=str, default=None, help='Name of cell group to be merged')
parser.add_argument('--chrom', type=str, default=None, help='Chromosome to impute')
parser.add_argument('--res', type=int, default=None, help='Bin size as integer to generate contact matrix')
parser.add_argument('--impute_mode', type=str, default=None, help='Suffix of imputed matrix file names')
parser.add_argument('--norm_mode', type=str, default='dist_trim', help='Suffix of normalized file names')

parser.add_argument('--min_dist', type=int, default=50000, help='Minimum distance threshold of loop')
parser.add_argument('--max_dist', type=int, default=10000000, help='Maximum distance threshold of loop')
parser.add_argument('--pad', type=int, default=5, help='One direction size of larger square for donut background')
parser.add_argument('--gap', type=int, default=2, help='One direction size of smaller square for donut background')
parser.add_argument('--thres_bl', type=int, default=1.33, help='Fold change threshold against bottom left background')
parser.add_argument('--thres_d', type=int, default=1.33, help='Fold change threshold against donut background')
parser.add_argument('--thres_h', type=int, default=1.2, help='Fold change threshold against horizontal background')
parser.add_argument('--thres_v', type=int, default=1.2, help='Fold change threshold against vertical background')
opt = parser.parse_args()

if opt.chrom[:3]=='chr':
    c = opt.chrom[3:]
else:
    c = opt.chrom

thres = stats.norm(0, 1).isf(0.025)
celllist = np.loadtxt(indir + 'merged/' + opt.cell_list, dtype=np.str)
tot = len(celllist)
Q = load_npz(opt.indir + 'chr' + c + '/' + opt.cell + '_chr' + c + '_' + opt.impute_mode + '.npz')
Qsum, Esum, Osum, Nsum, N2sum = [csr_matrix(Q.shape) for i in range(5)]
start_time = time.time()
for i,cell in enumerate(celllist):
	Q = load_npz(opt.indir + 'chr' + c + '/' + opt.cell + '_chr' + c + '_' + opt.impute_mode + '.npz')
	E = load_npz(opt.indir + 'chr' + c + '/' + opt.cell + '_chr' + c + '_' + opt.impute_mode + '_' + opt.norm_mode + '.E.npz')
	O = E.copy()
	O.data = (O.data > thres).astype(int)
	Qsum += Q
	Osum += O
	Esum += E

Qsum.data = Qsum.data / tot
Esum.data = Esum.data / tot
Osum.data = Osum.data / tot

save_npz(opt.indir + 'merged/' + opt.group + '_' + opt.impute_mode + '_' + opt.norm_mode + '.chr' + c + '.Q.npz', Qsum)
save_npz(opt.indir + 'merged/' + opt.group + '_' + opt.impute_mode + '_' + opt.norm_mode + '.chr' + c + '.E.npz', Esum)
save_npz(opt.indir + 'merged/' + opt.group + '_' + opt.impute_mode + '_' + opt.norm_mode + '.chr' + c + '.O.npz', Osum)
print('Merge cell', time.time() - start_time)

E = Esum.toarray()
O = Osum.toarray()
del Qsum, Esum, Osum

start_time = time.time()
oefilter = np.logical_and(E>0, O>0.1)
loop = np.where(oefilter)
distfilter = np.logical_and((loop[1] - loop[0]) > (opt.min_dist / opt.res), (loop[1] - loop[0]) < (opt.max_dist / opt.res))
loop = (loop[0][distfilter], loop[1][distfilter])

start_time = time.time()
eloop = np.zeros((tot, len(loop[0])))
for i,cell in enumerate(celllist):
	eloop[i] = load_npz(opt.indir + 'chr' + c + '/' + opt.cell + '_chr' + c + '_' + opt.impute_mode + '_' + opt.norm_mode + '.T.npz')[loop].A.ravel()

print('Load loop', time.time() - start_time)

pvr = np.array([stats.wilcoxon(xx, alternative='greater')[1] for xx in eloop.T])
pvt = stats.ttest_1samp(eloop, 0, axis=0)
pvt[1][pvt[0]>0] *= 2
pvt[1][pvt[0]<=0] = 1
pvt = pvt[1]
print('Test loop', time.time() - start_time)
del eloop

w = opt.pad * 2 + 1
start_time = time.time()

kernel_bl = np.zeros((w, w), np.float32)
kernel_bl[-opt.pad:, :(opt.pad - opt.gap)] = 1
kernel_bl[-(opt.pad - opt.gap):, :opt.pad] = 1

kernel_donut = np.ones((w, w), np.float32)
kernel_donut[opt.pad, :] = 0
kernel_donut[:, opt.pad] = 0
kernel_donut[(opt.pad - opt.gap):(opt.pad + opt.gap + 1), (opt.pad - opt.gap):(opt.pad + opt.gap + 1)] = 0

kernel_lr = np.ones((3, w), np.float32)
kernel_lr[:, (opt.pad - opt.gap):(opt.pad + opt.gap + 1)] = 0

kernel_bu = np.ones((w, 3), np.float32)
kernel_bu[(opt.pad - opt.gap):(opt.pad + opt.gap + 1), :] = 0

kernel_bl = kernel_bl / np.sum(kernel_bl)
kernel_donut = kernel_donut / np.sum(kernel_donut)
kernel_lr = kernel_lr / np.sum(kernel_lr)
kernel_bu = kernel_bu / np.sum(kernel_bu)

Ebl = cv2.filter2D(E, -1, kernel=kernel_bl) * (E>0)
Edonut = cv2.filter2D(E, -1, kernel=kernel_donut) * (E>0)
Elr = cv2.filter2D(E, -1, kernel=kernel_lr) * (E>0)
Ebu = cv2.filter2D(E, -1, kernel=kernel_bu) * (E>0)

bkfilter = np.logical_and(np.logical_and(E/Ebl > opt.thres_bl, E/Edonut > opt.thres_d), 
						np.logical_and(E/Elr > opt.thres_h, E/Ebu > opt.thres_v))
del Ebl, Edonut, Elr, Ebu, E

data = np.array([loop[0], loop[1], bkfilter[loop].astype(int), pvr, pvt]).T
np.save(opt.indir + 'merged/' + opt.group + '_' + opt.impute_mode + '_' + opt.norm_mode + '.chr' + c + '.loop.npy', data)
print('Filter loop', time.time() - start_time)

start_time = time.time()
Q = Qsum.tocoo()
data = np.array([np.zeros(len(Q.data)).astype(int), 
				np.repeat(['chr'+c], len(Q.data)), 
				Q.row * opt.res, 
				np.zeros(len(Q.data)).astype(int), 
				np.ones(len(Q.data)).astype(int), 
				np.repeat(['chr'+c], len(Q.data)), 
				Q.col * opt.res, 
				np.ones(len(Q.data)).astype(int), 
				np.around(Q.data * 100, decimals=4)]).T
data = pd.DataFrame(data, columns=['str1', 'chr1', 'x1', 'frag1','str2', 'chr2', 'y1', 'frag2', 'score'])
data.to_csv(opt.indir + 'merged/' + opt.group + '_' + opt.impute_mode + '_' + opt.norm_mode + '.chr' + c + '.Q.txt.gz', index=False, header=None, sep='\t', compression='gzip')
print('Write Q', time.time() - start_time)

start_time = time.time()
O = coo_matrix(O)
data = np.array([np.zeros(len(O.data)).astype(int), 
				np.repeat(['chr' + c], len(O.data)), 
				O.row * opt.res, 
				np.zeros(len(O.data)).astype(int), 
				np.ones(len(O.data)).astype(int), 
				np.repeat(['chr' + c], len(O.data)), 
				O.col * opt.res, 
				np.ones(len(O.data)).astype(int), 
				np.around(O.data * 100, decimals=2)]).T
data = pd.DataFrame(data, columns=['str1', 'chr1', 'x1', 'frag1','str2', 'chr2', 'y1', 'frag2', 'score'])
data.to_csv(opt.indir + 'merged/' + opt.group + '_' + opt.impute_mode + '_' + opt.norm_mode + '.chr' + c + '.O.txt.gz', index=False, header=None, sep='\t', compression='gzip')
print('Write O', time.time() - start_time)

'''
start_time = time.time()
E = coo_matrix(E)
data = np.array([np.zeros(len(E.data)).astype(int), np.repeat(['chr'+c], len(E.data)), E.row*res, np.zeros(len(E.data)).astype(int), np.ones(len(E.data)).astype(int), np.repeat(['chr'+c], len(E.data)), E.col*res, np.ones(len(E.data)).astype(int), np.around(E.data, decimals=2)]).T
data = pd.DataFrame(data, columns=['str1', 'chr1', 'x1', 'frag1','str2', 'chr2', 'y1', 'frag2', 'score'])
data.to_csv(indir + 'merged/' + ct + '_pad2_std1_rp0.5_sqrtvc_distnz_trim5.chr' + c + '.E.txt.gz', index=False, header=None, sep='\t', compression='gzip')
print('Write E', time.time() - start_time)
'''
