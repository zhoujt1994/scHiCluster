# command time python /gale/ddn/snm3C/humanPFC/code/concat_comp.py --indir /gale/ddn/snm3C/humanPFC/smoothed_matrix/${res0}b_resolution/ --cell_list /gale/ddn/snm3C/humanPFC/smoothed_matrix/celllist_long.txt --chrom ${SGE_TASK_ID} --mode raw --ncpus 10

import argparse
import numpy as np
from multiprocessing import Pool

def load_cpg(i):
        data = np.load(f'{opt.indir}chr{c}/{celllist[i]}_chr{c}_{opt.mode}.cpgcomp.npy')
        return [i, data[:-3], data[-3:]]

parser = argparse.ArgumentParser()
parser.add_argument('--indir', type=str, default=None, help='Directory of single-cell compartment score end with /')
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
result = p.map(load_cpg, np.arange(len(celllist)))
p.close()
comp = np.zeros((len(celllist), len(result[0][1])))
score = np.zeros((len(celllist), 3))
for i,x,s in result:
	comp[i] = x.copy()
	score[i] = s.copy()

np.save(f'{opt.indir}/merged/{opt.mode}_chr{c}.cpgcomp.npy', comp)
np.save(f'{opt.indir}/merged/{opt.mode}_chr{c}.cpgcompstr.npy', score)
