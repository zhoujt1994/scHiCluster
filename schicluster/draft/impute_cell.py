# command time python /gale/ddn/snm3C/humanPFC/code/impute_cell.py --indir /gale/raidix/rdx-5/zhoujt/projects/methylHiC/PFC_batch_merged/smoothed_matrix/1cell/${res0}b_resolution/chr${c}/ --outdir /gale/ddn/snm3C/humanPFC/smoothed_matrix/${res0}b_resolution/chr${c}/ --cell ${sample} --chrom ${c} --res ${res} --chrom_file /gale/netapp/home/zhoujt/genome/hg19/hg19.autosomal.chrom.sizes --mode pad2_std1_rp0.5_sqrtvc

import os
import time
import h5py
import cv2
cv2.useOptimized()
import argparse
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, save_npz, diags, eye
from scipy.sparse.linalg import norm

def impute_cell(indir, outdir, cell, chrom, res, chrom_file, 
                logscale=False, pad=1, std=1, rp=0.5, tol=0.01, 
                output_dist=500000000, output_format='hdf5', mode=None):

    def random_walk_cpu(P, rp, tol):
        if rp==1:
            return P
        I = eye(ngene)
        Q = P.copy()
        start_time = time.time()
        for i in range(30):
            Q_new = P.dot(Q * (1 - rp) + rp * I)
            delta = norm(Q - Q_new)
            Q = Q_new.copy()
            sparsity = Q.nnz / ngene / ngene
            end_time = time.time()
            print('Iter', i+1, 'takes', end_time-start_time, 'seconds, loss', delta, 'sparsity', sparsity)
            if delta < tol:
                break
        return Q

    if not os.path.exists(outdir):
        print('Output directory does not exist')
        return

    if chrom[:3]=='chr':
        c = chrom
    else:
        c = 'chr' + chrom
    if not mode:
        mode = f'pad{str(pad)}_std{str(std)}_rp{str(rp)}_sqrtvc'

    chromsize = pd.read_csv(chrom_file, sep='\t', header=None, index_col=0).to_dict()[1]

    start_time = time.time()
    ngene = int(chromsize[c] // res) + 1
    D = np.loadtxt(f'{indir}{cell}_{c}.txt')
    # to avoid bugs on chromosomes with 0/1 read
    if len(D)==0:
        D = np.array([[0,0,0]])
    elif len(D.shape)==1:
        D = D.reshape(1,-1)

    A = csr_matrix((D[:, 2], (D[:, 0], D[:, 1])), shape = (ngene, ngene))
    if logscale:
        A.data = np.log2(A.data + 1)

    end_time = time.time()
    print('Loading takes', end_time-start_time, 'seconds')

    start_time = time.time()
    B = cv2.GaussianBlur((A + A.T).astype(np.float32).toarray(), (pad*2+1, pad*2+1), std)
    end_time = time.time()
    print('Convolution takes', end_time-start_time, 'seconds')

    start_time = time.time()
    # remove diagonal before rwr
    B = csr_matrix(B)
    B = B - diags(B.diagonal())
    B = B + diags((B.sum(axis=0).A.ravel()==0).astype(int))
    d = diags(1 / B.sum(axis=0).A.ravel())
    P = d.dot(B)
    E = random_walk_cpu(P, rp, tol)
    print('RWR takes', time.time() - start_time, 'seconds')

    start_time = time.time()
    E = E + E.T
    d = E.sum(axis=0).A.ravel()
    d[d==0] = 1
    b = diags(1 / np.sqrt(d))
    E = b.dot(E).dot(b)
    print('SQRTVC takes', time.time() - start_time, 'seconds')

    # longest distance filter mask
    start_time = time.time()
    if (output_dist // res + 1) < ngene:
        idx = np.triu_indices(E.shape[0], 0)
        idxfilter = ((idx[1] - idx[0]) < (output_dist // res + 1))
        idx = (idx[0][idxfilter], idx[1][idxfilter])
        mask = csr_matrix((np.ones(len(idx[0])), (idx[0], idx[1])), E.shape)
        E = E.tocsr().multiply(mask)
    print('Filter takes', time.time() - start_time, 'seconds')

    if output_format=='npz':
        save_npz(f'{outdir}{cell}_{c}_{mode}.npz', E.astype(np.float32))
    else:
        f = h5py.File(f'{outdir}{cell}_{c}_{mode}.hdf5', 'w')
        g = f.create_group('Matrix')
        g.create_dataset('data', data=E.data, dtype='float32', compression='gzip')
        g.create_dataset('indices', data=E.indices, dtype=int, compression='gzip') 
        g.create_dataset('indptr', data=E.indptr, dtype=int, compression='gzip')
        g.attrs['shape'] = E.shape
        f.close()
    return

parser = argparse.ArgumentParser()
parser.add_argument('--indir', type=str, default=None, help='Directory of the contact matrix')
parser.add_argument('--outdir', type=str, default=None, help='Output directory end with /')
parser.add_argument('--cell', type=str, default=None, help='Specific identifier of a cell')
parser.add_argument('--chrom', type=str, default=None, help='Chromosome to impute')
parser.add_argument('--res', type=int, default=None, help='Bin size as integer to generate contact matrix')
parser.add_argument('--chrom_file', type=str, default=None, help='Path to the chromosome size files containing all chromosomes to be analyzed')

parser.add_argument('--logscale', dest='logscale', action='store_true', help='To log transform raw count')
parser.add_argument('--pad', type=int, default=1, help='Gaussian kernal size')
parser.add_argument('--std', type=float, default=1, help='Gaussian kernal standard deviation')

parser.add_argument('--rp', type=float, default=0.5, help='Restart probability of RWR')
parser.add_argument('--tol', type=float, default=0.01, help='Convergence tolerance of RWR')

parser.add_argument('--output_dist', type=int, default=500000000, help='Maximum distance threshold of contacts when writing output file')
parser.add_argument('--output_format', type=str, default='hdf5', help='Output file format (hdf5 or npz)')
parser.add_argument('--mode', type=str, default=None, help='Suffix of output file name')
opt = parser.parse_args()

impute_cell(opt.indir, opt.outdir, opt.cell, opt.chrom, opt.res, opt.chrom_file, 
            opt.logscale, opt.pad, opt.std, opt.rp, opt.tol, 
            opt.output_dist, opt.output_format, opt.mode)
