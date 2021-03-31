# command time python /gale/ddn/snm3C/humanPFC/code/loop_mergechr.py --inprefix /gale/ddn/snm3C/humanPFC/smoothed_matrix/10kb_resolution/merged/L23_pad1_std1_rp0.5_sqrtvc_dist_trim --outprefix /gale/ddn/snm3C/humanPFC/smoothed_matrix/10kb_resolution/merged/L23_pad1_std1_rp0.5_sqrtvc_dist_trim --chrom_file /gale/netapp/home/zhoujt/genome/hg19/hg19.autosomal.chrom.sizes 

import time
import argparse
import numpy as np
import pandas as pd
from heapq import *
from statsmodels.sandbox.stats.multicomp import multipletests as FDR

def loop_mergechr(inprefix, outprefix, chrom_file, split_file, res=10000, 
				thres_bl=1.33, thres_d=1.33, thres_h=1.2, thres_v=1.2, fdr_thres=0.1, dist_thres=20000, size_thres=1):

	def find_summit(loop, dist_thres):

		start_time = time.time()
		cord = loop[['x1', 'y1']].values // res
		idx = np.argsort(cord[:, 0])
		neighbor = {i:[] for i in range(len(idx))}
		for i in range(len(idx)-1):
			tmp = cord[idx[i]]
			for j in range(i+1,len(idx)):
				if cord[idx[j], 0] - tmp[0] > dist_thres:
					break
				if np.abs(tmp[1] - cord[idx[j], 1]) <= dist_thres:
					neighbor[idx[i]].append(idx[j])
					neighbor[idx[j]].append(idx[i])
		# dist = np.max([pairwise_distances(loop['x1'].values.reshape(-1,1)//res), 
		#				pairwise_distances(loop['y1'].values.reshape(-1,1)//res)], axis=0)
		# print('Compute distances takes', time.time() - start_time, 'seconds')
		# start_time = time.time()
		# G = nx.from_numpy_array((dist <= dist_thres))
		print('Build graph takes', time.time() - start_time, 'seconds')

		start_time = time.time()

		nodescore = loop['E'].values
		flag = np.zeros(len(nodescore))
		tot = len(nodescore)
		summit = []
		nodeheap = (loop['E'] * -1).reset_index().reset_index()[['E', 'level_0']].values.tolist()
		heapify(nodeheap)

		while tot>0:
			t = int(heappop(nodeheap)[1])
			while flag[t]:
				t = int(heappop(nodeheap)[1])
			q = [t]
			flag[t] = 1
			tot -= 1
			head = 0
			flagtmp = np.zeros(len(nodescore))
			while (head < len(q)):
				for t in neighbor[q[head]]:
					if not flagtmp[t] and nodescore[t]<nodescore[q[head]]:
						if not flag[t]:
							flag[t] = 1
							tot -= 1
						flagtmp[t] = 1
						q.append(t)
				head += 1
			summit.append([q[0], len(q)])
		summit = np.array(summit)
		loop = loop.iloc[summit[:,0]]
		loop['size'] = summit[:,1]
		print('BFS takes', time.time() - start_time, 'seconds')
		return loop

	chrom = np.loadtxt(chrom_file, dtype=np.str)[:,0]
	if not split_file:
		chrom_split = chrom.copy() 
	else:
		splitbed = pd.read_csv(split_file, sep='\t', header=None, index_col=0)
		chrom_split = np.concatenate([[c+'p', c+'q'] if c in splitbed.index else [c] for c in chrom])
	chrom_split = np.array([c if c[:3]=='chr' else 'chr'+c for c in chrom_split])

	start_time = time.time()
	loopall = []
	for c in chrom_split:
		data = pd.read_hdf(f'{inprefix}_{c}.loop.hdf5', key='loop')
		data['bkfilter'] = (((data['E']/data['E_bl'] > thres_bl) | (data['E_bl']<0)) & ((data['E']/data['E_donut'] > thres_d) | (data['E_donut']<0)) & ((data['E']/data['E_h'] > thres_h) | (data['E_h']<0)) & ((data['E']/data['E_v'] > thres_v) | (data['E_v']<0)))
		data['x1'] = data['x1'].astype(int) * res
		data['y1'] = data['y1'].astype(int) * res
		data['x2'] = data['x1'] + res
		data['y2'] = data['y1'] + res
		if c[-1]=='p' or c[-1]=='q':
			data['chr'] = c[:-1]
			if c[-1]=='q':
				data[['x1','x2','y1','y2']] += splitbed.loc[c,2] // res * res
		else:
			data['chr'] = c
		loopall.append(data)

	loopall = pd.concat(loopall, axis=0)
	loopall['rFDR'] = FDR(loopall['rpv'], 0.1, 'fdr_bh')[1]
	loopall['tFDR'] = FDR(loopall['tpv'], 0.1, 'fdr_bh')[1]
	print('Merge chromosomes takes', time.time() - start_time, 'seconds')

	loop = loopall.loc[(loopall['bkfilter']==1) & (loopall['tFDR']<fdr_thres)]
	loop.sort_values(by=['chr', 'x1', 'y1'])[['chr', 'x1', 'x2', 'chr', 'y1', 'y2', 'E']].to_csv(f'{outprefix}.loop.bedpe', sep='\t', index=False, header=None)
	scloop = loopall.loc[loopall['tFDR']<fdr_thres]
	scloop.sort_values(by=['chr', 'x1', 'y1'])[['chr', 'x1', 'x2', 'chr', 'y1', 'y2', 'E']].to_csv(f'{outprefix}.scloop.bedpe', sep='\t', index=False, header=None)
	bkloop = loopall.loc[loopall['bkfilter']==1]
	bkloop.sort_values(by=['chr', 'x1', 'y1'])[['chr', 'x1', 'x2', 'chr', 'y1', 'y2', 'E']].to_csv(f'{outprefix}.bkloop.bedpe', sep='\t', index=False, header=None)

	summit = pd.concat([find_summit(loop[loop['chr']==c], dist_thres//res) for c in chrom], axis=0)
	suumit = summit[summit['size'] >= size_thres]
	summit.sort_values(by=['chr', 'x1', 'y1'])[['chr', 'x1', 'x2', 'chr', 'y1', 'y2', 'E', 'size']].to_csv(f'{outprefix}.loopsummit.bedpe', sep='\t', index=False, header=None)

	return

parser = argparse.ArgumentParser()
parser.add_argument('--inprefix', type=str, default=None, help='Full path of a file containing the full path of all imputed files to be merged without .hdf5 suffix')
parser.add_argument('--outprefix', type=str, default=None, help='Prefix of merged matrix including directory')
parser.add_argument('--chrom_file', type=str, default=None, help='Path to the chromosome size files containing all chromosome to be analyzed as the first column') ## mm10 or hg38
parser.add_argument('--split_file', type=str, default=None, help='Path to the bed file containing all chromosomes need to split and one region per chromosome as splitting point')
parser.add_argument('--res', type=int, default=10000, help='Bin size as integer to generate contact matrix')
parser.add_argument('--thres_bl', type=int, default=1.33, help='Lowest fold change threshold against bottom left background')
parser.add_argument('--thres_d', type=int, default=1.33, help='Lowest fold change threshold against donut background')
parser.add_argument('--thres_h', type=int, default=1.2, help='Lowest fold change threshold against horizontal background')
parser.add_argument('--thres_v', type=int, default=1.2, help='Lowest fold change threshold against vertical background')
parser.add_argument('--fdr_thres', type=int, default=0.1, help='Highest t-test FDR threshold of loops')
parser.add_argument('--dist_thres', type=int, default=20000, help='Highest distance threshold to merge loops into summit')
parser.add_argument('--size_thres', type=int, default=1, help='Lowest loop number threshold of summit')
opt = parser.parse_args()

loop_mergechr(opt.inprefix, opt.outprefix, opt.chrom_file, opt.split_file, opt.res, 
		opt.thres_bl, opt.thres_d, opt.thres_h, opt.thres_v, opt.fdr_thres, opt.dist_thres, opt.size_thres)

