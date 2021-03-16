import time
import cv2
import argparse
import numpy as np
import pandas as pd
from scipy import stats
from scipy.sparse import save_npz, load_npz, csr_matrix, coo_matrix

parser = argparse.ArgumentParser()


def merge_cell(indir,
               cell_list,
               group,
               chrom,
               res,
               impute_mode,
               norm_mode='dist_trim',
               min_dist=50000,
               max_dist=10000000,
               pad=5,
               gap=2):

    if chrom[:3] == 'chr':
        c = chrom[3:]
    else:
        c = chrom

    thres = stats.norm(0, 1).isf(0.025)
    celllist = np.loadtxt(indir + 'merged/' + cell_list, dtype=np.str)
    tot = len(celllist)
    Q = load_npz(indir + 'chr' + c + '/' + celllist[0] + '_chr' + c + '_' + impute_mode + '.npz')
    Qsum, Esum, Osum, Nsum, N2sum = [csr_matrix(Q.shape) for i in range(5)]
    start_time = time.time()
    for i, cell in enumerate(celllist):
        Q = load_npz(indir + 'chr' + c + '/' + cell + '_chr' + c + '_' + impute_mode + '.npz')
        E = load_npz(
            indir + 'chr' + c + '/' + cell + '_chr' + c + '_' + impute_mode + '_' + norm_mode + '.E.npz')
        O = E.copy()
        O.data = (O.data > thres).astype(int)
        Qsum += Q
        Osum += O
        Esum += E

    Qsum.data = Qsum.data / tot
    Esum.data = Esum.data / tot
    Osum.data = Osum.data / tot

    save_npz(indir + 'merged/' + group + '_' + impute_mode + '_' + norm_mode + '.chr' + c + '.Q.npz',
             Qsum)
    save_npz(indir + 'merged/' + group + '_' + impute_mode + '_' + norm_mode + '.chr' + c + '.E.npz',
             Esum)
    save_npz(indir + 'merged/' + group + '_' + impute_mode + '_' + norm_mode + '.chr' + c + '.O.npz',
             Osum)
    print('Merge cell', time.time() - start_time)

    E = Esum.toarray()
    O = Osum.toarray()
    del Qsum, Esum, Osum

    start_time = time.time()
    oefilter = np.logical_and(E > 0, O > 0.1)
    loop = np.where(oefilter)
    distfilter = np.logical_and((loop[1] - loop[0]) > (min_dist / res),
                                (loop[1] - loop[0]) < (max_dist / res))
    loop = (loop[0][distfilter], loop[1][distfilter])

    start_time = time.time()
    eloop = np.zeros((tot, len(loop[0])))
    for i, cell in enumerate(celllist):
        eloop[i] = load_npz(
            indir + 'chr' + c + '/' + cell + '_chr' + c + '_' + impute_mode + '_' + norm_mode + '.T.npz')[
            loop].A.ravel()

    print('Load loop', time.time() - start_time)

    pvr = np.array([stats.wilcoxon(xx, alternative='greater')[1] for xx in eloop.T])
    pvt = stats.ttest_1samp(eloop, 0, axis=0)
    pvt[1][pvt[0] > 0] *= 2
    pvt[1][pvt[0] <= 0] = 1
    pvt = pvt[1]
    print('Test loop', time.time() - start_time)
    del eloop

    w = pad * 2 + 1
    start_time = time.time()

    kernel_bl = np.zeros((w, w), np.float32)
    kernel_bl[-pad:, :(pad - gap)] = 1
    kernel_bl[-(pad - gap):, :pad] = 1

    kernel_donut = np.ones((w, w), np.float32)
    kernel_donut[pad, :] = 0
    kernel_donut[:, pad] = 0
    kernel_donut[(pad - gap):(pad + gap + 1), (pad - gap):(pad + gap + 1)] = 0

    kernel_lr = np.ones((3, w), np.float32)
    kernel_lr[:, (pad - gap):(pad + gap + 1)] = 0

    kernel_bu = np.ones((w, 3), np.float32)
    kernel_bu[(pad - gap):(pad + gap + 1), :] = 0

    kernel_bl = kernel_bl / np.sum(kernel_bl)
    kernel_donut = kernel_donut / np.sum(kernel_donut)
    kernel_lr = kernel_lr / np.sum(kernel_lr)
    kernel_bu = kernel_bu / np.sum(kernel_bu)

    Ebl = cv2.filter2D(E, -1, kernel=kernel_bl) * (E > 0)
    Edonut = cv2.filter2D(E, -1, kernel=kernel_donut) * (E > 0)
    Elr = cv2.filter2D(E, -1, kernel=kernel_lr) * (E > 0)
    Ebu = cv2.filter2D(E, -1, kernel=kernel_bu) * (E > 0)

    data = np.array([loop[0], loop[1], pvr, pvt, (E/Ebl)[loop], (E/Edonut)[loop], (E/Elr)[loop], (E/Ebu)[loop]]).T
    np.save(indir + 'merged/' + group + '_' + impute_mode + '_' + norm_mode + '.chr' + c + '.loop.npy',
            data)
    print('Filter loop', time.time() - start_time)

    return


'''
start_time = time.time()
E = coo_matrix(E)
data = np.array([np.zeros(len(E.data)).astype(int), np.repeat(['chr'+c], len(E.data)), 
E.row*res, np.zeros(len(E.data)).astype(int), np.ones(len(E.data)).astype(int), 
np.repeat(['chr'+c], len(E.data)), E.col*res, np.ones(len(E.data)).astype(int), np.around(E.data, decimals=2)]).T
data = pd.DataFrame(data, columns=['str1', 'chr1', 'x1', 'frag1','str2', 'chr2', 'y1', 'frag2', 'score'])
data.to_csv(indir + 'merged/' + ct + '_pad2_std1_rp0.5_sqrtvc_distnz_trim5.chr' + c + '.E.txt.gz', 
index=False, header=None, sep='\t', compression='gzip')
print('Write E', time.time() - start_time)
'''
