import os
import sys
import time
import numpy as np
import torch
import torch.nn.functional as F
from scipy.sparse import csr_matrix
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans

def neighbor_ave_gpu(A, pad):
	if pad==0:
		return torch.from_numpy(A).float().cuda()
	ll = pad * 2 + 1
	conv_filter = torch.ones(1, 1, ll, ll).cuda()
	B = F.conv2d(torch.from_numpy(A[None, None,: ,:]).float().cuda(), conv_filter, padding = pad * 2)
	return (B[0, 0 ,pad:-pad, pad:-pad] / float(ll * ll))

def random_walk_gpu(A, rp):
	ngene, _ = A.shape
	A = A - torch.diag(torch.diag(A))
	A = A + torch.diag(torch.sum(A, 0) == 0).float()
	P = torch.div(A, torch.sum(A, 0))
	Q = torch.eye(ngene).cuda()
	I = torch.eye(ngene).cuda()
	for i in range(30):
		Q_new = (1 - rp) * I + rp * torch.mm(Q, P)
		delta = torch.norm(Q - Q_new, 2)
		Q = Q_new
		if delta < 1e-6:
			break
	return Q

def neighbor_ave_cpu(A, pad):
	if pad==0:
		return A
	ngene, _ = A.shape
	ll = pad * 2 + 1
	B, C, D, E = [np.zeros((ngene + ll, ngene + ll)) for i in range(4)]
	B[(pad + 1):(pad + ngene + 1), (pad + 1):(pad + ngene + 1)] = A[:]
	F = B.cumsum(axis = 0).cumsum(axis = 1)
	C[ll :, ll:] = F[:-ll, :-ll]
	D[ll:, :] = F[:-ll, :]
	E[:, ll:] = F[:, :-ll]
	return (np.around(F + C - D - E, decimals=8)[ll:, ll:] / float(ll * ll))

def random_walk_cpu(A, rp):
	ngene, _ = A.shape
	A = A - np.diag(np.diag(A))
	A = A + np.diag(np.sum(A, axis=0) == 0)
	P = np.divide(A, np.sum(A, axis = 0))
	Q = np.eye(ngene)
	I = np.eye(ngene)
	for i in range(30):
		Q_new = (1 - rp) * I + rp * np.dot(Q, P)
		delta = np.linalg.norm(Q - Q_new)
		Q = Q_new.copy()
		if delta < 1e-6:
			break
	return Q

def hicluster_gpu(network, chromsize, nc, res=1000000, pad=1, rp=0.5, prct=20, ndim=20):
	matrix=[]
	for i, c in enumerate(chromsize):
		ngene = int(chromsize[c] / res)+1
		start_time = time.time()
		Q_concat = torch.zeros(len(label), ngene * ngene).float().cuda()
		j = 0
		for cell in network:
			D = np.loadtxt(cell + '_chr' + c + '.txt')
			A = csr_matrix((D[:, 2], (D[:, 0], D[:, 1])), shape = (ngene, ngene)).toarray()
			A = np.log2(A + A.T + 1)
			A = neighbor_ave_gpu(A, pad)
			Q = random_walk_gpu(A, rp)
			Q_concat[j, :] = ((Q > np.percentile(Q, 100 - prct)) * (Q < 1.0)).reshape(ngene * ngene)
			j += 1
		end_time = time.time()
		print('Load and impute chromosome', c, 'take', end_time - start_time, 'seconds')
		ndim = int(min(Q_concat.shape) * 0.2) - 1
		# U, S, V = torch.svd(Q_concat, some=True)
		# R_reduce = torch.mm(U[:, :ndim], torch.diag(S[:ndim])).cuda().numpy()
		Q_concat = Q_concat.cpu().numpy()
		pca = PCA(n_components = ndim)
		R_reduce = pca.fit_transform(Q_concat)
		matrix.append(R_reduce)
		print(c)
	matrix = np.concatenate(matrix, axis = 1)
	pca = PCA(n_components = min(matrix.shape) - 1)
	matrix_reduce = pca.fit_transform(matrix)
	kmeans = KMeans(n_clusters = nc, n_init = 200).fit(matrix_reduce[:, :ndim])
	return kmeans.labels_, matrix_reduce

def hicluster_cpu(network, chromsize, nc, res=1000000, pad=1, rp=0.5, prct=20, ndim=20):
	matrix=[]
	for i, c in enumerate(chromsize):
		ngene = int(chromsize[c] / res)+1
		start_time = time.time()
		Q_concat = np.zeros((len(network), ngene * ngene)).astype(float)
		j = 0
		for cell in network:
			D = np.loadtxt(cell + '_chr' + c + '.txt')
			A = csr_matrix((D[:, 2], (D[:, 0], D[:, 1])), shape = (ngene, ngene)).toarray()
			A = np.log2(A + A.T + 1)
			A = neighbor_ave_cpu(A, pad)
			Q = random_walk_cpu(A, rp)
			Q_concat[j, :] = ((Q > np.percentile(Q, 100 - prct)) * (Q < 1.0)).reshape(ngene * ngene)
			j += 1
		end_time = time.time()
		print('Load and impute chromosome', c, 'take', end_time - start_time, 'seconds')
		ndim = int(min(Q_concat.shape) * 0.2) - 1
		pca = PCA(n_components = ndim)
		R_reduce = pca.fit_transform(Q_concat)
		matrix.append(R_reduce)
		print(c)
	matrix = np.concatenate(matrix, axis = 1)
	pca = PCA(n_components = min(matrix.shape) - 1)
	matrix_reduce = pca.fit_transform(matrix)
	kmeans = KMeans(n_clusters = nc, n_init = 200).fit(matrix_reduce[:, :ndim])
	return kmeans.labels_, matrix_reduce

def pca(network, chromsize, nc, res=1000000, ndim=20):
	matrix=[]
	for i, c in enumerate(chromsize):
		ngene = int(chromsize[c] / res)+1
		start_time = time.time()
		uptri = np.triu_indices(ngene, 1)
		A_concat = np.zeros((len(label), len(uptri[0]))).astype(float)
		j = 0
		for cell in network:
			D = np.loadtxt(cell + '_chr' + c + '.txt')
			A = csr_matrix((D[:, 2], (D[:, 0], D[:, 1])), shape = (ngene, ngene)).toarray()
			A = np.log2(A + A.T + 1)
			A_concat[j, :] = A[uptri]
			j += 1
		end_time = time.time()
		print('Load and impute chromosome', c, 'take', end_time - start_time, 'seconds')
		pca = PCA(n_components = min(A_concat.shape)-1)
		A_reduce = pca.fit_transform(A_concat)
		matrix.append(A_reduce)
		print(c)
	matrix = np.concatenate(matrix, axis = 1)
	pca = PCA(n_components = min(matrix.shape) - 1)
	matrix_reduce = pca.fit_transform(matrix)
	kmeans = KMeans(n_clusters = nc, n_init = 200).fit(matrix_reduce[:, 1:(ndim + 1)])
	return kmeans.labels_, matrix_reduce

def compartment(network, chromsize, nc, res=1000000, ndim=20):
	global chrcg
	pca = PCA(n_components = 2)
	matrix = []
	for i, c in enumerate(chrom):
		ngene = int(chromsize[c] / res)+1
		start_time = time.time()
		comp = np.zeros((len(label), ngene)).astype(float)
		j = 0
		for cell in network:
			D = np.loadtxt(cell+'_chr'+c+'.txt')
			A = csr_matrix((D[:, 2], (D[:, 0], D[:, 1])), shape=(ngene, ngene)).toarray()
			A = np.log2(A + A.T + 1)
			B = np.corrcoef(A)
			B[np.isnan(B)] = 0
			pc = zscore(pca.fit_transform(B)[:,0])
			if np.mean(chrcg[c][pc>0])>np.mean(chrcg[c][pc<0]):
				pc = -pc.copy()
			comp[j] = pc.copy()
			j += 1
		end_time = time.time()
		print('Load and random walk take', end_time - start_time, 'seconds')
		matrix.append(comp)
		print(c)
	matrix = np.concatenate(matrix, axis = 1)
	pca = PCA(n_components = min(matrix.shape) - 1)
	matrix_reduce = pca.fit_transform(matrix)
	kmeans = KMeans(n_clusters = nc, n_init = 200).fit(matrix_reduce[:, :ndim])
	return kmeans.labels_, matrix_reduce



cg = np.loadtxt('/cellar/users/zhoujt1994/genome/hg19/bin/hg19.1mb.bin.CpG.txt', dtype=np.str, skiprows=1, usecols=(0,9,11,12))
cgdata = cg[:,1:].astype(float)
cgdata = cgdata[:,2]/(cgdata[:,1]-cgdata[:,0])
cgdata[np.isnan(cgdata)] = 0.0
chrcg = {c:cgdata[cg[:,0]=='chr'+c] for c in chrom}
ctlist = ['HeLa', 'HAP1', 'GM12878', 'K562']
network = {ct:np.loadtxt('/cellar/users/zhoujt1994/projects/scHiC/Ramani2017/cell_matrix/1mb_resolution/'+ct+'/sample_list.txt', dtype=np.str) for ct in ctlist}
label = np.concatenate([[ct for i in range(len(network[ct]))] for ct in ctlist])
network = np.concatenate([network[ct] for ct in ctlist])
mm9dim = [197195432,181748087,159599783,155630120,152537259,149517037,152524553,131738871,124076172,129993255,121843856,121257530,120284312,125194864,103494974,98319150,95272651,90772031,61342430]
hg19dim = [249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566,155270560]
chrom = [str(i+1) for i in range(22)] + ['X']
chromsize = {chrom[i]:hg19dim[i] for i in range(len(chrom))}
start_time = time.time()
cluster, embedding = compartment(network, chromsize, nc=4)
print(time.time() - start_time)
print(ARI(label, cluster))


56.43739366531372
75.48667979240417
39.40321183204651
51.541046380996704
381065.965

cg = np.loadtxt('/cellar/users/zhoujt1994/genome/mm9/bin/mm9.1mb.bin.CpG.txt', dtype=np.str, skiprows=1, usecols=(0,9,11,12))
cgdata = cg[:,1:].astype(float)
cgdata = cgdata[:,2]/(cgdata[:,1]-cgdata[:,0])
cgdata[np.isnan(cgdata)] = 0.0
chrcg = {c:cgdata[cg[:,0]=='chr'+c] for c in chrom}
ctlist = ['NSN', 'SN', 'ZygP', 'ZygM']
network = {ct:np.loadtxt('/cellar/users/zhoujt1994/projects/scHiC/Flyamer2017/cell_matrix/1mb_resolution/'+ct+'/sample_list.txt', dtype=np.str) for ct in ctlist}
label = np.concatenate([[ct for i in range(len(network[ct]))] for ct in ctlist]).astype('U8')
label[label == 'NSN'] = 'Oocyte'
label[label == 'SN'] = 'Oocyte'
network = np.concatenate([network[ct] for ct in ctlist])
mm9dim = [197195432,181748087,159599783,155630120,152537259,149517037,152524553,131738871,124076172,129993255,121843856,121257530,120284312,125194864,103494974,98319150,95272651,90772031,61342430]
chrom = [str(i+1) for i in range(19)]
chromsize = {chrom[i]:mm9dim[i] for i in range(len(chrom))}
start_time = time.time()
cluster, embedding = pca(network, chromsize, nc=3)
print(time.time() - start_time)
kmeans = KMeans(n_clusters = 3, n_init = 200).fit(embedding[:, :20])
print(ARI(label, kmeans.labels_))


25.54651141166687
30.359964847564697
20.048206090927124
29.159992933273315
31783.053

