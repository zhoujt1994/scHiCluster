# command time python /gale/ddn/snm3C/humanPFC/code/concat_ins.py --indir /gale/ddn/snm3C/humanPFC/smoothed_matrix/${res0}b_resolution/ --cell_list /gale/ddn/snm3C/humanPFC/smoothed_matrix/celllist_long.txt --chrom ${SGE_TASK_ID} --mode pad1_std1_rp0.5_sqrtvc --ncpus 10

import argparse
import numpy as np
from multiprocessing import Pool

def load_ins(i):
	data = np.load(f'{opt.indir}chr{c}/{celllist[i]}_chr{c}_{opt.mode}.is.npy')
	return [i, data]

parser = argparse.ArgumentParser()
parser.add_argument('--indir', type=str, default=None, help='Directory of single-cell insulation score end with /')
parser.add_argument('--cell_list', type=str, default=None, help='Full path of a file containing a list of cell identifiers to be concatenate')
parser.add_argument('--chrom', type=str, default=None, help='Chromosome to impute')
parser.add_argument('--mode', type=str, default=None, help='Suffix of imputed matrix file names')
parser.add_argument('--ncpus', type=int, default=10, help='# threads for parallelization')
opt = parser.parse_args()

if opt.chrom[:3]=='chr':
	c = opt.chrom[3:]
else:
	c = opt.chrom
celllist = np.loadtxt(opt.cell_list, dtype=np.str)
p = Pool(opt.ncpus)
result = p.map(load_ins, np.arange(len(celllist)))
p.close()
ins = np.zeros((len(celllist), len(result[0][1])))
for i,x in result:
	ins[i] = x.copy()

np.save(f'{opt.indir}/merged/{opt.mode}_chr{c}.ins.npy', ins)
