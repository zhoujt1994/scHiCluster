# command time python /gale/ddn/snm3C/humanPFC/code/concat_ins.py --cell_list /gale/ddn/snm3C/humanPFC/smoothed_matrix/filelist/inslist_pad1_std1_rp0.5_sqrtvc_chr${c}.txt --outprefix /gale/ddn/snm3C/humanPFC/smoothed_matrix/${res0}b_resolution/merged/pad1_std1_rp0.5_sqrtvc_chr${c} --ncpus 10

import argparse
import numpy as np
from multiprocessing import Pool

def load_ins(i):
	data = np.load(celllist[i])
	return [i, data]

parser = argparse.ArgumentParser()
parser.add_argument('--cell_list', type=str, default=None, help='Full path of a file containing the full path of all insulation npy files to be concatenate')
parser.add_argument('--outprefix', type=str, default=None, help='Prefix of concatenated matrix including directory')
parser.add_argument('--ncpus', type=int, default=10, help='# threads for parallelization')
opt = parser.parse_args()

celllist = np.loadtxt(opt.cell_list, dtype=np.str)
p = Pool(opt.ncpus)
result = p.map(load_ins, np.arange(len(celllist)))
p.close()
ins = np.zeros((len(celllist), len(result[0][1])))
for i,x in result:
	ins[i] = x.copy()

np.save(f'{opt.outprefix}.ins.npy', ins)
