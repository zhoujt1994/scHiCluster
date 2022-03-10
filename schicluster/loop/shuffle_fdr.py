import cooler
import numpy as np
import pandas as pd
from scipy.sparse import load_npz, csr_matrix, save_npz, triu
from scipy.stats import rankdata


def _t_score(T, T2, tot):
    t = T2 - T.power(2)
    t.data = np.sqrt((tot - 1) / t.data)
    t = T.multiply(t)
    t.data[np.isnan(t.data)] = 0
    return t


def compute_t(group_prefix, tot=None):
    cool_t = cooler.Cooler(f'{group_prefix}.E.cool')
    cool_t2 = cooler.Cooler(f'{group_prefix}.E2.cool')
    if not tot:
        tot = cool_t.info['group_n_cells']
    for chrom in cool_t.chromnames:
        T = triu(cool_t.matrix(balance=False, sparse=True).fetch(chrom))
        T2 = triu(cool_t2.matrix(balance=False, sparse=True).fetch(chrom))
        t = _t_score(T, T2, tot)
        save_npz(f'{group_prefix}_{chrom}.tglobal.npz', t)

    cool_t = cooler.Cooler(f'{group_prefix}.T.cool')
    cool_t2 = cooler.Cooler(f'{group_prefix}.T2.cool')
    for chrom in cool_t.chromnames:
        T = triu(cool_t.matrix(balance=False, sparse=True).fetch(chrom))
        T2 = triu(cool_t2.matrix(balance=False, sparse=True).fetch(chrom))
        t = _t_score(T, T2, tot)
        save_npz(f'{group_prefix}_{chrom}.tlocal.npz', t)
    return tot


def permute_fdr(chrom_size_path,
                black_list_path,
                shuffle_group_prefix,
                real_group_prefix,
                res=10000,
                pad=7,
                min_dist=5,
                max_dist=500):
    chrom_size_series = pd.read_csv(chrom_size_path,
                                    sep='\t',
                                    index_col=0,
                                    header=None,
                                    squeeze=True)
    chrom_size = (chrom_size_series.values // res).astype(int) + 1
    chroms = chrom_size_series.index
    covfilter = []
    cumsize = [[0] for _ in range(min_dist, max_dist)]
    bkl = pd.read_csv(black_list_path, sep='\t', header=None, index_col=None)
    for k, chrom in enumerate(chroms):
        cov = np.zeros(chrom_size[k])
        for xx, yy in bkl.loc[bkl[0] == chrom, [1, 2]].values // res:
            cov[max([xx - pad, 0]):min([len(cov), yy + pad + 1])] = 1
        tmp = []
        row, col = np.diag_indices(chrom_size[k])
        for i in range(min_dist, max_dist):
            tmp.append(np.logical_or(cov[row[:-i]], cov[col[i:]]))
            cumsize[i - min_dist].append(cumsize[i - min_dist][-1] +
                                         np.sum(~tmp[-1]))
        covfilter.append(tmp)

    for bktype in ['local', 'global']:
        t = [[] for _ in range(min_dist, max_dist)]
        tnull = [[] for _ in range(min_dist, max_dist)]
        for k, chrom in enumerate(chroms):
            tmp = load_npz(
                f'{shuffle_group_prefix}_{chrom}.t{bktype}.npz').toarray()
            for i in range(min_dist, max_dist):
                tmp0 = tmp.diagonal(i).copy()[~covfilter[k][i - min_dist]]
                tnull[i - min_dist].append(tmp0)
            tmp = load_npz(
                f'{real_group_prefix}_{chrom}.t{bktype}.npz').toarray()
            for i in range(min_dist, max_dist):
                tmp0 = tmp.diagonal(i).copy()[~covfilter[k][i - min_dist]]
                t[i - min_dist].append(tmp0)

        fdr = []
        for i in range(len(t)):
            tmp1 = -np.concatenate(t[i])
            tmp2 = -np.concatenate(tnull[i])
            t1 = rankdata(np.concatenate((tmp1, tmp2)),
                          method='max')[:len(tmp1)]
            t2 = rankdata(tmp1, method='max')
            fdr.append((t1 - t2) / t2)

        for k, chrom in enumerate(chroms):
            ngene = chrom_size[k]
            tmp = np.zeros((ngene, ngene))
            row, col = np.diag_indices(ngene)
            for i in range(min_dist, max_dist):
                idx = (row[:-i][~covfilter[k][i - min_dist]],
                       col[i:][~covfilter[k][i - min_dist]])
                tmp[row[:-i], col[i:]] = 1
                tmp[idx] = fdr[i - min_dist][(cumsize[i - min_dist][k]):(
                    cumsize[i - min_dist][k + 1])]
            tmp = csr_matrix(tmp)
            save_npz(f'{shuffle_group_prefix}_{chrom}.permutefdr{bktype}.npz',
                     tmp)
    return


def update_fdr_qval(chrom_size_path,
                    real_group_prefix,
                    shuffle_group_prefix,
                    res=10000,
                    min_dist=5,
                    max_dist=500):
    chrom_size_series = pd.read_csv(chrom_size_path,
                                    sep='\t',
                                    index_col=0,
                                    header=None,
                                    squeeze=True)
    data: pd.DataFrame = pd.read_hdf(f'{real_group_prefix}.totalloop_info.hdf',
                                     key='data')
    data['global_qval'] = 1
    data['local_qval'] = 1
    data = data.loc[((data['distance'] // res) > min_dist)
                    & ((data['distance'] // res) < max_dist)
                    & data['bkfilter']]
    for chrom in chrom_size_series.index:
        tmpfilter = (data['chrom'] == chrom)
        tmp = data.loc[tmpfilter, ['x1', 'y1']].values // res
        coord = (tmp[:, 0], tmp[:, 1])
        tmp = load_npz(f'{shuffle_group_prefix}_{chrom}.permutefdrlocal.npz')
        data.loc[tmpfilter, 'local_qval'] = tmp[coord].A.ravel()
        tmp = load_npz(f'{shuffle_group_prefix}_{chrom}.permutefdrglobal.npz')
        data.loc[tmpfilter, 'global_qval'] = tmp[coord].A.ravel()

    data.to_hdf(f'{real_group_prefix}.totalloop_info.hdf',
                key='data',
                format='table')
    return data
