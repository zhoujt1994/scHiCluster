import time
import cooler
import numpy as np
from scipy.ndimage import convolve
from scipy.sparse import csr_matrix, save_npz, triu
from scipy.stats import zscore


def calc_diag_stats(E, n_dims):
    """Calculate cutoff, average, std, count of non-zero pixels of each diagonals of the E"""
    ave, std, top, count = np.zeros((4, n_dims), dtype=np.float32)
    for i in range(n_dims):
        tmp = E.diagonal(i)
        if tmp.size == 0:
            top[i] = 0
            ave[i] = 0
            std[i] = 0
            count[i] = 0
        else:
            cutoff = np.percentile(tmp, 99)
            tmp = np.where(tmp < cutoff, tmp, cutoff)
            top[i] = cutoff
            ave[i] = np.mean(tmp)
            std[i] = np.std(tmp)
            count[i] = np.sum(tmp > 0)
        # TODO smoothing
    return ave, std, top, count


def calculate_chrom_background_normalization(cell_url,
                                             chrom,
                                             resolution,
                                             output_prefix,
                                             dist=10050000,
                                             cap=5,
                                             pad=5,
                                             gap=2,
                                             min_cutoff=1e-6,
                                             log_e=False):
    """
    Compute the background for each chromosome in each cell

    Parameters
    ----------
    cell_url
    chrom
    resolution
    output_prefix
    dist
    cap
    pad
    gap
    min_cutoff
    log_e

    Returns
    -------
    E is the global diagonal normalized matrix
    T is the local background normalized version of E
    """
    cell_cool = cooler.Cooler(cell_url)
    # Load the cell imputed matrix as E
    E = triu(cell_cool.matrix(balance=False, sparse=True).fetch(chrom))
    E = E.astype(np.float32).toarray()

    # create an upper triangle mask
    mask = np.zeros(E.shape, dtype=bool)
    row, col = np.diag_indices(E.shape[0])
    mask[row, col] = True
    for i in range(1, dist // resolution + 1):
        mask[row[:-i], col[i:]] = True

    if log_e:
        # normalize E at log scale
        E[row, col] = 0
        for i in range(1, dist // resolution + 1):
            tmp = E.diagonal(i).copy()
            tmp_filter = (tmp > 0)
            tmp2 = tmp[tmp_filter]
            if len(tmp2) == 0:
                E[row[:-i], col[i:]] = 0
            else:
                tmp2 = zscore(np.log10(tmp2))
                tmp[~tmp_filter] = tmp2.min()
                tmp[tmp_filter] = tmp2.copy()
                E[row[:-i], col[i:]] = tmp.copy()
    else:
        # normalize E at linear scale
        # Calculate the diagonal stats of E
        n_dims = dist // resolution + pad + 1
        ave, std, top, count = calc_diag_stats(E, n_dims)
        # normalize E with the diagonal backgrounds
        idx = np.where(mask)
        tmp = E[idx]
        tmp = np.min([tmp, top[idx[1] - idx[0]]], axis=0)
        tmp = (tmp - ave[idx[1] - idx[0]]) / (std[idx[1] - idx[0]] + 1e-5)  # add a small value to prevent divide by 0
        tmp[count[idx[1] - idx[0]] < 100] = 0
        tmp[std[idx[1] - idx[0]] == 0] = 0
        tmp[tmp > cap] = cap
        tmp[tmp < -cap] = -cap
        E[idx] = tmp.copy()
    # print(f'Norm {time.time() - start_time:.3f}', E.dtype, tmp.dtype)

    # normalize E with the local backgrounds to generate T
    w = pad * 2 + 1
    kernel = np.ones((w, w), np.float32)
    kernel[(pad - gap):(pad + gap + 1), (pad - gap):(pad + gap + 1)] = 0
    kernel = kernel / np.sum(kernel)
    T = convolve(E, kernel, mode='mirror')
    E = csr_matrix(E)
    T = csr_matrix(T * mask)
    if min_cutoff > 0:
        # mask out small abs values
        E = E.multiply(np.abs(E) > min_cutoff)
        T = T.multiply(np.abs(T) > min_cutoff)
    T = E - T
    # print(f'Bkg {time.time() - start_time:.3f}', E.dtype, T.dtype)

    save_npz(f'{output_prefix}.E.npz', E)
    save_npz(f'{output_prefix}.T.npz', T)
    return
