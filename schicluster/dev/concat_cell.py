import time
import h5py
import numpy as np
from scipy.sparse import save_npz, csr_matrix, vstack
from sklearn.decomposition import TruncatedSVD


def run_svd(matrix):
    svd = TruncatedSVD(n_components=50, algorithm='arpack')
    matrix_reduce = svd.fit_transform(matrix)
    matrix_reduce = matrix_reduce / svd.singular_values_
    return matrix_reduce


def concat_cell(indir, chrom, mode, res, dist=10000000, save_raw=True, scale_factor=100000):
    # TODO simplify chrom definition
    if chrom[:3] == 'chr':
        chrom = chrom[3:]

    # TODO simplify file IO here
    celllist = np.loadtxt(indir + '../cell_list.txt', dtype=np.str)

    f = h5py.File(indir + 'chr' + chrom + '/' + celllist[0] + '_chr' + chrom + '_' + mode + '.hdf5', 'r')
    ngene = f['Matrix'].attrs['shape'][0]
    f.close()
    idx = np.triu_indices(ngene, k=1)
    idxfilter = np.array([(yy - xx) < (dist / res + 1) for xx, yy in zip(idx[0], idx[1])])
    idx = (idx[0][idxfilter], idx[1][idxfilter])

    start_time = time.time()
    matrix = []
    for i, cell in enumerate(celllist):
        f = h5py.File(indir + 'chr' + chrom + '/' + cell + '_chr' + chrom + '_' + mode + '.hdf5', 'r')
        g = f['Matrix']
        A = csr_matrix((g['data'][()], g['indices'][()], f['indptr'][()]), g.attrs['shape'])
        f.close()
        matrix.append(csr_matrix(A[idx]))
        if i % 100 == 0:
            print(i, 'cells loaded', time.time() - start_time, 'seconds')

    matrix = vstack(matrix)

    if save_raw:
        save_npz(indir + 'merged/' + mode + '_chr' + chrom + '.npz', matrix)

    matrix.data = matrix.data * scale_factor
    matrix_reduce = run_svd(matrix.data)
    np.save(indir + 'merged/' + mode + '_chr' + chrom + '.svd50.npy', matrix_reduce)
    return
