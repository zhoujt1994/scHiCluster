# command time python /gale/ddn/snm3C/humanPFC/code/loop_sumcell_chr.py --cell_list /gale/ddn/snm3C/humanPFC/smoothed_matrix/10kb_resolution/filelist/L23_covgroup${i}_pad1_std1_rp0.5_sqrtvc_chr${c}_looplist.txt --outprefix /gale/ddn/snm3C/humanPFC/smoothed_matrix/10kb_resolution/merged/L23_covgroup${i}_pad1_std1_rp0.5_sqrtvc_dist_trim_chr${c} --res 10000
# command time python /gale/ddn/snm3C/humanPFC/code/loop_sumcell_chr.py --cell_list /gale/ddn/snm3C/humanPFC/smoothed_matrix/10kb_resolution/filelist/L23_pad1_std1_rp0.5_sqrtvc_chr${c}_looplist.txt --outprefix /gale/ddn/snm3C/humanPFC/smoothed_matrix/10kb_resolution/merged/L23_pad1_std1_rp0.5_sqrtvc_dist_trim_chr${c} --res 10000 --grouplist /gale/ddn/snm3C/humanPFC/smoothed_matrix/10kb_resolution/filelist/L23_pad1_std1_rp0.5_sqrtvc_chr${c}_grouplist.txt

import time
import cv2
import h5py
import argparse
import numpy as np
import pandas as pd
from scipy import stats
from scipy.sparse import save_npz, load_npz, csr_matrix, coo_matrix

def loop_sumcell_chr(cell_list, outprefix, res, group_list=None, matrix='QEO',
            norm_mode='dist_trim', min_dist=50000, max_dist=10000000, pad=5, gap=2):

    thres = stats.norm(0, 1).isf(0.025)
    celllist = np.loadtxt(cell_list, dtype=np.str)
    tot = len(celllist)
    with h5py.File(f'{celllist[0]}.hdf5', 'r') as f:
        dim1, dim2 = f['Matrix'].attrs['shape']
    Qsum, Esum, Osum = [csr_matrix((dim1, dim2)) for i in range(3)]
    start_time = time.time()
    if group_list:
        grouplist = pd.read_csv(group_list, sep='\t', index_col=0, header=None)
        for i, cell in enumerate(grouplist.index):
            if 'Q' in matrix:
                Q = load_npz(f'{cell}.Q.npz')
                Q.data = Q.data * grouplist.loc[cell][1]
                Qsum += Q
            if 'E' in matrix:
                E = load_npz(f'{cell}.E.npz')
                E.data = E.data * grouplist.loc[cell][1]
                Esum += E
            if 'O' in matrix:
                O = load_npz(f'{cell}.O.npz')
                O.data = O.data * grouplist.loc[cell][1]
                Osum += O
            print('Merge', i, 'groups takes', time.time() - start_time, 'seconds')
    else:
        for i, cell in enumerate(celllist):
            if 'Q' in matrix:
                with h5py.File(f'{cell}.hdf5', 'r') as f:
                    g = f['Matrix']
                    Q = csr_matrix((g['data'][()], g['indices'][()], g['indptr'][()]), g.attrs['shape'])
                Qsum += Q
            if ('E' in matrix) or ('O' in matrix):
                E = load_npz(f'{cell}_{norm_mode}.E.npz')
            if 'E' in matrix:
                Esum += E
            if 'O' in matrix:
                O = E.copy()
                O.data = (O.data > thres).astype(int)
                Osum += O
            print('Merge', i, 'cells takes', time.time() - start_time, 'seconds')
    if 'Q' in matrix:
        Qsum.data = Qsum.data / tot
        save_npz(f'{outprefix}.Q.npz', Qsum)
    if 'E' in matrix:
        Esum.data = Esum.data / tot
        save_npz(f'{outprefix}.E.npz', Esum)
    if 'O' in matrix:
        Osum.data = Osum.data / tot
        save_npz(f'{outprefix}.O.npz', Osum)
    print('Merge cell', time.time() - start_time)

    if (not 'E' in matrix) or (not 'O' in matrix):
        return

    E = Esum.toarray()
    O = Osum.toarray()
    del Qsum, Esum, Osum

    start_time = time.time()
    oefilter = np.logical_and(E > 0, O > 0.1)
    loop = np.where(oefilter)
    distfilter = np.logical_and((loop[1] - loop[0]) > (min_dist / res), (loop[1] - loop[0]) < (max_dist / res))
    loop = (loop[0][distfilter], loop[1][distfilter])

    start_time = time.time()
    eloop = np.zeros((len(celllist), len(loop[0])))
    for i, cell in enumerate(celllist):
        eloop[i] = load_npz(f'{cell}_{norm_mode}.T.npz')[loop].A.ravel()
    print('Load', len(loop[0]), 'loop', time.time() - start_time)

    start_time = time.time()
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

    data = pd.DataFrame(np.array([loop[0], loop[1], pvr, pvt, E[loop], Ebl[loop], Edonut[loop], Elr[loop], Ebu[loop]]).T, columns = ['x1', 'y1', 'rpv', 'tpv', 'E', 'E_bl', 'E_donut', 'E_h', 'E_v'])
    data.to_hdf(f'{outprefix}.loop.hdf5', key='loop')
    print('Bulk bkg', time.time() - start_time)

    return
'''
parser = argparse.ArgumentParser()
parser.add_argument('--cell_list', type=str, default=None, help='Full path of a file containing the full path of all imputed files to be merged without .hdf5 suffix')
parser.add_argument('--outprefix', type=str, default=None, help='Prefix of merged matrix including directory')
parser.add_argument('--res', type=int, default=None, help='Bin size as integer to generate contact matrix')
parser.add_argument('--group_list', type=str, default=None, help='Full path of a file containing the full path of all merged files to be merged')
parser.add_argument('--matrix', type=str, default='QEO', help='Types of matrices to merge') #Q E O
parser.add_argument('--norm_mode', type=str, default='dist_trim', help='Suffix of normalized file names')
parser.add_argument('--min_dist', type=int, default=50000, help='Minimum distance threshold of loop')
parser.add_argument('--max_dist', type=int, default=10000000, help='Maximum distance threshold of loop')
parser.add_argument('--pad', type=int, default=5, help='One direction size of larger square for donut background')
parser.add_argument('--gap', type=int, default=2, help='One direction size of smaller square for donut background')
opt = parser.parse_args()

loop_sumcell_chr(opt.cell_list, opt.outprefix, opt.res, opt.group_list, opt.matrix,
        opt.norm_mode, opt.min_dist, opt.max_dist, opt.pad, opt.gap)
'''
