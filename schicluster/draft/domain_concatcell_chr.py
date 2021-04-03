# for c in `seq 1 22`; do awk -v c=$c '{printf("/gale/ddn/snm3C/humanPFC/smoothed_matrix/25kb_resolution/chr%s/%s_chr%s_pad2_std1_rp0.5_sqrtvc.w10.domain.bed\n",c,$1,c)}' celllist_long.txt > 25kb_resolution/filelist/domainlist_pad2_std1_rp0.5_sqrtvc_chr${c}.txt; echo $c; done
# for c in `seq 1 22`; do awk -v c=$c '{printf("/gale/ddn/snm3C/humanPFC/smoothed_matrix/25kb_resolution/chr%s/%s_chr%s_pad2_std1_rp0.5_sqrtvc.w10.ins.npy\n",c,$1,c)}' celllist_long.txt > 25kb_resolution/filelist/inslist_pad2_std1_rp0.5_sqrtvc_chr${c}.txt; echo $c; done

# command time python /gale/ddn/snm3C/humanPFC/code/domain_concatcell_chr.py --cell_list /gale/ddn/snm3C/humanPFC/smoothed_matrix/25kb_resolution/filelist/inslist_pad2_std1_rp0.5_sqrtvc_chr${c}.txt --outprefix /gale/ddn/snm3C/humanPFC/smoothed_matrix/${res0}b_resolution/merged/pad2_std1_rp0.5_sqrtvc_chr${c}.w10 --res 25000 --input_type insulation --ncpus 10
# command time python /gale/ddn/snm3C/humanPFC/code/domain_concatcell_chr.py --cell_list /gale/ddn/snm3C/humanPFC/smoothed_matrix/25kb_resolution/filelist/domainlist_pad2_std1_rp0.5_sqrtvc_chr${c}.txt --outprefix /gale/ddn/snm3C/humanPFC/smoothed_matrix/${res0}b_resolution/merged/pad2_std1_rp0.5_sqrtvc_chr${c}.w10 --res 25000 --input_type boundary --ncpus 10

import argparse
import numpy as np
import pandas as pd
from multiprocessing import Pool
from scipy.sparse import csr_matrix, save_npz

def load_insulation(i):
	try:
		data = np.load(celllist[i])
	except:
		print(celllist[i])
		data = 0
	return [i, data]

def load_boundary(i):
	try:
		tmp = pd.read_csv(celllist[i], sep='\t', header=None)
	except:
		print(celllist[i])
		data = 0
	else:
		data = np.zeros(int(tmp.iloc[-1,2] // rs) + 1)
		tmp = tmp[tmp[3]=='domain'][[1,2]].values // rs
		data[tmp[:,0]] += 1
		data[tmp[:,1]] += 1
	return [i, data]

def domain_concatcell_chr(cell_list, outprefix, res, input_type='insulation', ncpus=10):

	global celllist, rs
	celllist = np.loadtxt(cell_list, dtype=np.str)
	rs = res

	p = Pool(ncpus)
	if input_type=='insulation':
		result = p.map(load_insulation, np.arange(len(celllist)))
	elif input_type=='boundary':
		result = p.map(load_boundary, np.arange(len(celllist)))
	ins = np.zeros((len(celllist), len(result[0][1])))
	for i,x in result:
		if not isinstance(x, int):
			ins[i] = x.copy()
	if input_type=='insulation':
		np.save(f'{outprefix}.{input_type}.npy', ins)
	elif input_type=='boundary':
		save_npz(f'{outprefix}.{input_type}.npz', csr_matrix(ins))
	p.close()
	return

parser = argparse.ArgumentParser()
parser.add_argument('--cell_list', type=str, default=None, help='Full path of a file containing the full path of all insulation npy or domain txt files to be concatenate')
parser.add_argument('--outprefix', type=str, default=None, help='Prefix of concatenated matrix including directory')
parser.add_argument('--res', type=int, default=None, help='Bin size as integer')
parser.add_argument('--input_type', type=str, default='insulation', help='Whether input files are insulation.npy or domain.txt') # insulation or boundary
parser.add_argument('--ncpus', type=int, default=10, help='# threads for parallelization')
opt = parser.parse_args()

domain_concatcell_chr(opt.cell_list, opt.outprefix, opt.res, opt.input_type, opt.ncpus)
