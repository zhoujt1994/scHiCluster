import time
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, diags, eye, save_npz
from scipy.sparse.linalg import norm
from scipy.ndimage import gaussian_filter
import cooler
import logging

# from ..cool import write_coo


def calc_sparsity(matrix):
    row, col = matrix.shape
    sparsity = matrix.nnz / row / col
    return sparsity


def random_walk_cpu(P, rp, tol):
    if rp == 1:
        return P

    _start_time = time.time()
    n_genes = P.shape[0]
    I = eye(n_genes, dtype=np.float32)
    Q = P.copy()
    for i in range(30):
        Q_new = P.dot(Q * (1 - rp) + rp * I)
        delta = norm(Q - Q_new)
        Q = Q_new.copy()
        sparsity = calc_sparsity(Q)
        _end_time = time.time()
        logging.debug(
            f'Iter {i + 1} takes {(_end_time - _start_time):.3f} seconds. '
            f'Loss: {delta:.3f}; Sparsity: {sparsity:.3f}', P.dtype, Q.dtype)
        if delta < tol:
            break
    return Q


def impute_chromosome(chrom,
                      resolution,
                      output_path,
                      scool_url=None,
                      contact_path=None,
                      chrom_size_path=None,
                      logscale=False,
                      pad=1,
                      std=1,
                      rp=0.5,
                      tol=0.01,
                      window_size=500000000,
                      step_size=10000000,
                      output_dist=500000000,
                      min_cutoff=0,
                      chrom1=1,
                      pos1=2,
                      chrom2=5,
                      pos2=6):
    if scool_url is not None:
        cell_cool = cooler.Cooler(scool_url)
        A = cell_cool.matrix(balance=False, sparse=True).fetch(chrom)
        # A = A + diags(A.diagonal())
        n_bins = A.shape[0]
    elif contact_path is not None:
        if chrom_size_path is not None:
            chrom_sizes = pd.read_csv(chrom_size_path, sep='\t', index_col=0, header=None).squeeze(axis=1)
            if chrom=='all':
                from schicluster.cool.utilities import get_chrom_offsets
                bins_df = cooler.binnify(chrom_sizes, resolution)
                chrom_offset = get_chrom_offsets(bins_df)
                n_bins = bins_df.shape[0]
            else:
                n_bins = (chrom_sizes.loc[chrom] // resolution) + 1
        else:
            print("ERROR : Must provide chrom_size_path if using contact file as input")
            return
        A = pd.read_csv(contact_path, sep='\t', header=None, index_col=None, comment='#')[[chrom1, pos1, chrom2, pos2]]
        if chrom=='all':
            A = A[A[chrom1].isin(chrom_offset) & A[chrom2].isin(chrom_offset)]
            A[pos1] = A[chrom1].map(chrom_offset) + (A[pos1] - 1) // resolution
            A[pos2] = A[chrom2].map(chrom_offset) + (A[pos2] - 1) // resolution
        else:
            A = A.loc[(A[chrom1]==chrom) & (A[chrom2]==chrom)]
            A[[pos1, pos2]] = (A[[pos1, pos2]] - 1) // resolution
        A = A.groupby(by=[pos1, pos2])[chrom1].count().reset_index()
        A = csr_matrix((A[chrom1].astype(np.int32), (A[pos1], A[pos2])), (n_bins, n_bins))
        A = A + A.T
    else:
        print("ERROR : Must provide either scool_url or contact_file_path")
        return

    ws = int(window_size // resolution)
    ss = int(step_size // resolution)

    # log transform
    if logscale:
        A.data = np.log2(A.data + 1)

    # Remove diagonal before convolution
    A = A - diags(A.diagonal())

    # Gaussian convolution and
    start_time = time.time()
    if pad > 0:
        # full matrix step
        A = gaussian_filter(A.astype(np.float32).toarray(),
                            std, order=0, mode='mirror', truncate=pad)
        A = csr_matrix(A)
    # else:
    #     A = A + A.T
    end_time = time.time()
    logging.debug(f'Convolution takes {end_time - start_time:.3f} seconds')

    # Remove diagonal before RWR
    A = A - diags(A.diagonal())

    # Random Walk with Restart
    start_time = time.time()
    if ws >= n_bins or rp == 1:
        B = A + diags((A.sum(axis=0).A.ravel() == 0).astype(int))
        d = diags(1 / B.sum(axis=0).A.ravel())
        P = d.dot(B).astype(np.float32)
        E = random_walk_cpu(P, rp, tol)
    else:
        # if the chromosome is too large, compute by chunks
        idx = (np.repeat(np.arange(ws), ws), np.tile(np.arange(ws), ws))
        idxfilter = (np.abs(idx[1] - idx[0]) < (output_dist // resolution + 1))
        idx = (idx[0][idxfilter], idx[1][idxfilter])
        # first filter
        idxfilter = ((idx[0] + idx[1]) < (ws + ss))
        idx1 = (idx[0][idxfilter], idx[1][idxfilter])
        mask1 = csr_matrix((np.ones(len(idx1[0])), (idx1[0], idx1[1])),
                           (ws, ws))
        # last filter
        idxfilter = ((idx[0] + idx[1]) >= (
                (n_bins - ws) // ss * 2 + 1) * ss + 3 * ws - 2 * n_bins)
        idx2 = (idx[0][idxfilter], idx[1][idxfilter])
        mask2 = csr_matrix((np.ones(len(idx2[0])), (idx2[0], idx2[1])),
                           (ws, ws))
        # center filter
        idxfilter = np.logical_and((idx[0] + idx[1]) < (ws + ss),
                                   (idx[0] + idx[1]) >= (ws - ss))
        idx0 = (idx[0][idxfilter], idx[1][idxfilter])
        mask0 = csr_matrix((np.ones(len(idx0[0])), (idx0[0], idx0[1])),
                           (ws, ws))

        start_time = time.time()
        E = csr_matrix(A.shape, dtype=np.float32)
        for ll in [x for x in range(0, n_bins - ws, ss)] + [n_bins - ws]:
            B = A[ll:(ll + ws), ll:(ll + ws)]
            B = B + diags((B.sum(axis=0).A.ravel() == 0).astype(int))
            d = diags(1 / B.sum(axis=0).A.ravel())
            P = d.dot(B).astype(np.float32)
            Etmp = random_walk_cpu(P, rp, tol)
            if ll == 0:
                E[ll:(ll + ws), ll:(ll + ws)] += Etmp.multiply(mask1)
            elif ll == (n_bins - ws):
                E[ll:(ll + ws), ll:(ll + ws)] += Etmp.multiply(mask2)
            else:
                E[ll:(ll + ws), ll:(ll + ws)] += Etmp.multiply(mask0)
    logging.debug(f'RWR takes {time.time() - start_time:.3f} seconds')

    # Normalize
    start_time = time.time()
    E += E.T
    d = E.sum(axis=0).A.ravel()
    d[d == 0] = 1
    b = diags(1 / np.sqrt(d))
    E = b.dot(E).dot(b)
    logging.debug(f'SQRTVC takes {time.time() - start_time:.3f} seconds')

    start_time = time.time()
    # mask the lower triangle of E
    # TODO This part is MEM intensive, the mask below can be combined with the chunk mask above
    idx = np.triu_indices(E.shape[0], 0)
    if (output_dist // resolution + 1) < n_bins:
        # longest distance filter mask
        idxfilter = ((idx[1] - idx[0]) < (output_dist // resolution + 1))
        idx = (idx[0][idxfilter], idx[1][idxfilter])
    mask = csr_matrix((np.ones(len(idx[0])), (idx[0], idx[1])),
                      E.shape,
                      dtype=np.float32)
    E = E.tocsr().multiply(mask)
    logging.debug(f'Filter takes {time.time() - start_time:.3f} seconds')

    # TODO put this part inside RWR, before normalize
    # min_cutoff = tol/
    # Make values < min_cutoff to 0
    if min_cutoff > 0:
        s_before = calc_sparsity(E)
        E = E.multiply(E > min_cutoff)
        s_after = calc_sparsity(E)
        logging.debug(f'Mask values smaller than {min_cutoff}. Sparsity before {s_before:.3f}, after {s_after:.3f}')

    # save to file
    # write_coo(output_path, E)
    save_npz(output_path, E)

    return
