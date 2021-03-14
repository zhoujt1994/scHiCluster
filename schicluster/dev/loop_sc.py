import time
import cv2

cv2.useOptimized()
import numpy as np
from scipy.sparse import load_npz, save_npz, csr_matrix


def loop_sc(outdir, cell, chrom, impute_mode, res, dist, cap, pad, gap, norm_mode):
    if chrom[:3] == 'chr':
        c = chrom[3:]
    else:
        c = chrom

    E = load_npz(outdir + cell + '_chr' + c + '_' + impute_mode + '.npz')
    start_time = time.time()
    ave, std, top, count = np.zeros((4, dist // res + pad + 1))
    for i in range(dist // res + pad + 1):
        tmp = E.diagonal(i)
        top[i] = np.percentile(tmp, 99)
        tmp[tmp > top[i]] = top[i]
        ave[i] = np.mean(tmp)
        std[i] = np.std(tmp)
        count[i] = np.sum(tmp > 0)

    print('Curve', time.time() - start_time, '#Nonzero', np.sum(count))
    start_time = time.time()
    E.data = np.min([E.data, top[E.col - E.row]], axis=0)
    E = E.astype(np.float32).toarray()

    idx = np.triu_indices(E.shape[0], 0)
    idxfilter = ((idx[1] - idx[0]) < (dist // res + 1))
    idx = (idx[0][idxfilter], idx[1][idxfilter])
    mask = csr_matrix((np.ones(len(idx[0])), (idx[0], idx[1])), E.shape)

    tmp = E[idx]
    tmp = (tmp - ave[idx[1] - idx[0]]) / std[idx[1] - idx[0]]
    tmp[count[idx[1] - idx[0]] < 100] = 0
    tmp[std[idx[1] - idx[0]] == 0] = 0
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

    save_npz(outdir + cell + '_chr' + c + '_' + impute_mode + '_' + norm_mode + '.E.npz', E)
    save_npz(outdir + cell + '_chr' + c + '_' + impute_mode + '_' + norm_mode + '.T.npz', N)
    return
