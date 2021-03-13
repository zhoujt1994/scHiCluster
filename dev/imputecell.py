import time
import h5py
import cv2
cv2.useOptimized()
import argparse
import numpy as np
from scipy.sparse import csr_matrix, save_npz, diags, eye
from scipy.sparse.linalg import norm

def random_walk_cpu(P, rp, tol, dist, spr):
    global ngene
    if rp==1:
        return P
    I = eye(ngene)
    Q = rp * I + (1 - rp) * P
    if dist < ngene:
        idx = np.triu_indices(ngene, 0)
        idxfilter = ((idx[1] - idx[0]) < dist)
        idx = (np.concatenate((idx[0][idxfilter], idx[1][idxfilter])), np.concatenate((idx[1][idxfilter], idx[0][idxfilter])))
        mask = csr_matrix((np.ones(len(idx[0])), (idx[0], idx[1])), P.shape)
    start_time = time.time()
    for i in range(30):
        Q_new = rp * I + (1 - rp) * P.dot(Q)
        delta = norm(Q - Q_new)
        Q = Q_new.copy()
        sparsity = Q.nnz / ngene / ngene
        if (dist < ngene) and (sparsity > spr):
            Q = Q.multiply(mask)
        end_time = time.time()
        print('Iter', i+1, 'takes', end_time-start_time, 'seconds, loss', delta, 'sparsity', sparsity)
        if delta < tol:
            break
    return Q

parser = argparse.ArgumentParser()
parser.add_argument('--indir', type=str, default=None, help='Directory of the contact matrix')
parser.add_argument('--outdir', type=str, default=None, help='Output directory end with /')
parser.add_argument('--cell', type=str, default=None, help='Specific identifier of a cell')
parser.add_argument('--chrom', type=str, default=None, help='Chromosome to impute')
parser.add_argument('--res', type=int, default=None, help='Bin size as integer to generate contact matrix')
parser.add_argument('--genome', type=str, default=None, help='Genome assembly version') ## mm10 or hg38

parser.add_argument('--logscale', type=bool, default=False, help='Whether to log transform raw count or not')
parser.add_argument('--pad', type=int, default=1, help='Gaussian kernal size')
parser.add_argument('--std', type=float, default=1, help='Gaussian kernal standard deviation')

parser.add_argument('--rp', type=float, default=0.5, help='Restart probability of RWR')
parser.add_argument('--tol', type=float, default=0.01, help='Convergence tolerance of RWR')
parser.add_argument('--rwr_dist', type=int, default=500000000, help='Maximum distance threshold of contacts when doing RWR')
parser.add_argument('--rwr_sparsity', type=float, default=1, help='Minimum sparsity to apply rwr_dist')

parser.add_argument('--output_dist', type=int, default=500000000, help='Maximum distance threshold of contacts when writing output file')
parser.add_argument('--output_format', type=str, default='hdf', help='Output file format (hdf5 or npz)')
parser.add_argument('--mode', type=str, default=None, help='Suffix of output file name')
opt = parser.parse_args()

if opt.chrom[:3]=='chr':
    c = opt.chrom[3:]
else:
    c = opt.chrom
if not opt.mode:
    opt.mode = 'pad' + str(opt.pad) + '_std' + str(opt.std) + '_rp' + str(opt.rp) + '_sqrtvc'

hg38dim = [248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717, 133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468]
hg19dim = [249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566]
mm10dim = [195471971, 182113224, 160039680, 156508116, 151834684, 149736546, 145441459, 129401213, 124595110, 130694993, 122082543, 120129022, 120421639, 124902244, 104043685, 98207768, 94987271, 90702639, 61431566]

if opt.genome=='hg38':
    chrom = [str(i+1) for i in range(22)]
    chromsize = {x:y for x,y in zip(chrom, hg38dim)}
elif opt.genome=='hg19':
    chrom = [str(i+1) for i in range(22)]
    chromsize = {x:y for x,y in zip(chrom, hg19dim)}
elif opt.genome=='mm10':
    chrom = [str(i+1) for i in range(19)]
    chromsize = {x:y for x,y in zip(chrom, mm10dim)}

start_time = time.time()
ngene = int(chromsize[c] // opt.res) + 1
D = np.loadtxt(opt.indir + opt.cell + '_chr' + c + '.txt')
if len(D)==0:
    D = np.array([[0,0,0]])
elif len(D.shape)==1:
    D = D.reshape(1,-1)

A = csr_matrix((D[:, 2], (D[:, 0], D[:, 1])), shape = (ngene, ngene))
if opt.logscale:
    A.data = np.log2(A.data + 1)

end_time = time.time()
print('Loading chr', c, 'takes', end_time-start_time, 'seconds')

start_time = time.time()
B = cv2.GaussianBlur((A + A.T).toarray().astype(np.float32), (opt.pad*2+1, opt.pad*2+1), opt.std)
end_time = time.time()
print('Convolution chr', c, 'take', end_time-start_time, 'seconds')

start_time = time.time()
B = csr_matrix(B)
B = B - diags(B.diagonal())
B = B + diags((B.sum(axis=0).A.ravel()==0).astype(int))
d = diags(1 / B.sum(axis=0).A.ravel())
P = d.dot(B)
Q = random_walk_cpu(P, opt.rp, opt.tol, int(opt.rwr_dist // opt.res) + 1, opt.rwr_sparsity)
print('RWR', time.time() - start_time)

start_time = time.time()
E = Q + Q.T
d = E.sum(axis=0).A.ravel()
d[d==0] = 1
b = diags(1 / np.sqrt(d))
E = b.dot(E).dot(b)
print('SQRTVC', time.time() - start_time)

start_time = time.time()
idx = np.triu_indices(E.shape[0], 0)
idxfilter = ((idx[1] - idx[0]) < (opt.output_dist // opt.res + 1))
idx = (idx[0][idxfilter], idx[1][idxfilter])
mask = csr_matrix((np.ones(len(idx[0])), (idx[0], idx[1])), E.shape)
E = E.multiply(mask)
E = E.tocsr()
print('Filter', time.time() - start_time)

if opt.output_format=='npz':
    save_npz(opt.outdir + opt.cell + '_chr' + c + '_' + opt.mode + '.npz', E)
else:
    f = h5py.File(opt.outdir + opt.cell + '_chr' + c + '_' + opt.mode + '.hdf5', 'w')
    g = f.create_group('Matrix')
    g.create_dataset('data', data=E.data, dtype='float32', compression='gzip')
    g.create_dataset('indices', data=E.indices, dtype=int, compression='gzip') 
    g.create_dataset('indptr', data=E.indptr, dtype=int, compression='gzip')
    g.attrs['shape'] = E.shape
    f.close()
