import os
import sys
import time
import numpy as np
import torch
import torch.nn.functional as F
from multiprocessing import Pool
from scipy.sparse import csr_matrix
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans, SpectralClustering
from sklearn.metrics.cluster import adjusted_rand_score as ARI

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

def impute_gpu(args):
	cell, c, ngene, pad, rp, prct = args
	D = np.loadtxt(cell + '_chr' + c + '.txt')
	A = csr_matrix((D[:, 2], (D[:, 0], D[:, 1])), shape = (ngene, ngene)).toarray()
	A = np.log2(A + A.T + 1)
	A = neighbor_ave_gpu(A, pad)
	if rp==-1:
		Q = A[:]
	else:
		Q = random_walk_gpu(A, rp)
	if prct==-1:
		Q = Q.reshape(ngene * ngene)
	else:
		Q = ((Q > np.percentile(Q, 100 - prct)) * (Q < 1.0)).reshape(ngene * ngene)
	return Q

def hicluster_gpu(network, chromsize, nc, res=1000000, pad=1, rp=0.5, prct=20, ndim=20):
	matrix=[]
	for i, c in enumerate(chromsize):
		ngene = int(chromsize[c] / res)+1
		start_time = time.time()
		Q_concat = torch.zeros(len(network), ngene * ngene).float().cuda()
		for j, cell in enumerate(network):
			Q_concat[j] = impute_gpu([cell, c, ngene, pad, rp, prct])
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

def impute_cpu(args):
	cell, c, ngene, pad, rp, prct = args
	D = np.loadtxt(cell + '_chr' + c + '.txt')
	A = csr_matrix((D[:, 2], (D[:, 0], D[:, 1])), shape = (ngene, ngene)).toarray()
	A = np.log2(A + A.T + 1)
	A = neighbor_ave_cpu(A, pad)
	if rp==-1:
		Q = A[:]
	else:
		Q = random_walk_cpu(A, rp)
	if prct==-1:
		Q = Q.reshape(ngene * ngene)
	else:
		Q = ((Q > np.percentile(Q, 100 - prct)) * (Q < 1.0)).reshape(ngene * ngene)
	return [cell, Q]

def hicluster_cpu(network, chromsize, nc, res=1000000, pad=1, rp=0.5, prct=20, ndim=20, ncpus=10):
	matrix=[]
	for i, c in enumerate(chromsize):
		ngene = int(chromsize[c] / res)+1
		start_time = time.time()
		paras = [[cell, c, ngene, pad, rp, prct] for cell in network]
		p = Pool(ncpus)
		result = p.map(impute, paras)
		p.close()
		index = {x[0]:j for j,x in enumerate(result)}
		Q_concat = np.array([result[index[x]][1] for x in network])
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

def raw_pca(network, chromsize, nc, res=1000000, ndim=20):
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

def ds_pca(network, chromsize, nc, res=1000000, ndim=20):
	matrix=[]
	for i, c in enumerate(chromsize):
		ngene = int(chromsize[c] / res)+1
		start_time = time.time()
		uptri = np.triu_indices(ngene, 1)
		A_concat = np.zeros((len(label), len(uptri[0]))).astype(float)
		j = 0
		tot = int(np.min([np.sum(np.loadtxt(cell + '_chr' + c + '.txt')[:, 2]) for cell in network]))
		for cell in network:
			D = np.loadtxt(cell + '_chr' + c + '.txt')		
			A = csr_matrix((ngene, ngene)).toarray().astype(float)
			idx = np.concatenate([[j for i in range(int(D[j,2]))] for j in range(len(D))])
			shu = np.arange(len(idx))
			np.random.shuffle(shu)
			for x in idx[shu[:tot]]:
				A[int(D[x,0]), int(D[x,1])] += 1
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
	kmeans = KMeans(n_clusters = nc, n_init = 200).fit(matrix_reduce[:, :ndim])
	return kmeans.labels_, matrix_reduce

def comp_comp(O, cg):
	ngene, _ = Q.shape
	N = torch.zeros(ngene, ngene).float().cuda()
	for k in range(ngene):
		diagsum = torch.sum(torch.diag(O, k)) / (ngene - k)
		N = N + torch.diag(torch.diag(O, k) / diagsum, k)
	D = torch.diag(torch.diag(N))
	N = N + torch.t(N) - D
	N[N!=N] = 0.0
	C = corrcoef(N).cpu().numpy()
	pca = PCA(n_components=2)
	pc = pca.fit_transform(C)
	corr = np.abs([np.corrcoef(cg, pc[:,i])[0,1] for i in range(2)])
	# if corr[0] > corr[1]:
	# 	k = 0
	# else:
	# 	k = 1
	k = 0
	pc, corr = pc[:,k], corr[k]
	comp = torch.mean(N[pc>0][:,pc<0]).item() * 2 / (torch.mean(N[pc>0][:,pc>0]).item() + torch.mean(N[pc<0][:,pc<0]).item())
	if np.corrcoef(cg, pc)[0,1] < 0:
		return [C, -pc, corr, comp]
	else:
		return [C, pc, corr, comp]

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
			A = neighbor_ave_gpu(A, 0)
			comp[j] = comp_comp(A, chrcg[c])
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

def decay(network, chromsize, nc, res=1000000, ndim=20):
	matrix = []
	for i, c in enumerate(chromsize):
		ngene = int(chromsize[c] / res)+1
		start_time = time.time()
		dec = np.zeros((len(label), ngene-1)).astype(float)
		j = 0
		for j, cell in enumerate(network):
			D = np.loadtxt(cell+'_chr'+c+'.txt')
			A = csr_matrix((D[:, 2], (D[:, 0], D[:, 1])), shape=(ngene, ngene)).toarray()
			tmp = np.array([np.sum(np.diag(A, k)) for k in range(1,ngene)])
			dec[j] = tmp / np.sum(tmp)
		end_time = time.time()
		print('Load and random walk take', end_time - start_time, 'seconds')
		matrix.append(dec)
		print(c)
	matrix = np.concatenate(matrix, axis = 1)
	pca = PCA(n_components = min(matrix.shape) - 1)
	matrix_reduce = pca.fit_transform(matrix)
	kmeans = KMeans(n_clusters = nc, n_init = 200).fit(matrix_reduce[:, :ndim])
	return kmeans.labels_, matrix_reduce

def merge_gpu(network, c, res, pad=1, rp=0.5, prct=-1):
	global chromsize
	ngene = int(chromsize[c] / res) + 1
	start_time = time.time()
	Q_sum = torch.zeros(ngene * ngene).float().cuda()
	for cell in network:
		Q_sum = Q_sum + impute_gpu([cell, c, ngene, pad, rp, prct])
	end_time = time.time()
	print('Load and impute chromosome', c, 'take', end_time - start_time, 'seconds')
	Q_sum = Q_sum.cpu().numpy().reshape(ngene, ngene)
	return Q_sum

def merge_cpu(network, c, res, pad=1, rp=0.5, prct=-1):
	global chromsize
	ngene = int(chromsize[c] / res) + 1
	start_time = time.time()
	Q_sum = torch.zeros(ngene * ngene).float().cuda()
	for cell in network:
		Q_sum = Q_sum + impute_cpu([cell, c, ngene, pad, rp, prct])[1]
	end_time = time.time()
	print('Load and impute chromosome', c, 'take', end_time - start_time, 'seconds')
	Q_sum = Q_sum.reshape(ngene, ngene)
	return Q_sum

def output_topdom(cell, c, Q, res):
	global chromsize
	ngene, _ = Q.shape
	B = [['chr' + c, i * res, (i + 1) * res] for i in range(ngene)]
	B[-1][-1] = chromsize[c]
	C = np.concatenate((B, Q), axis=1)
	np.savetxt(cell + '_chr' + c + '.topdommatrix', C, fmt = '%s', delimiter = '\t')
	return

def output_sparse(cell, c, Q, res):
	idx = np.where(Q>0)
	C = np.concatenate(((idx[0]*res)[:, None], (idx[1]*res)[:,None], Q[idx][:, None]), axis=1)
	np.savetxt(cell + '_chr' + c + '.sparsematrix', C, fmt = '%d\t%d\t%s', delimiter = '\n')
	return

mm9dim = [197195432,181748087,159599783,155630120,152537259,149517037,152524553,131738871,124076172,129993255,121843856,121257530,120284312,125194864,103494974,98319150,95272651,90772031,61342430]
hg19dim = [249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566,155270560]

# File list and labels of dataset Ramani 2017
ctlist = ['HeLa', 'HAP1', 'GM12878', 'K562']
network = [np.loadtxt('1mb_resolution/'+ct+'/sample_list.txt',dtype=np.str) for ct in ctlist]
label = np.array([ctlist[i] for i in range(len(ctlist)) for j in range(len(network[i]))]).astype('U8')
network = np.concatenate(network)
chrom = [str(i+1) for i in range(22)] + ['X']
chromsize = {chrom[i]:hg19dim[i] for i in range(len(chrom))}
nc = 4

# CpG content for each bin
cg = np.loadtxt('hg19/bin/hg19.1mb.bin.CpG.txt', dtype=np.str, skiprows=1, usecols=(0,9,11,12))
cgdata = cg[:,1:].astype(float)
cgdata = cgdata[:,2]/(cgdata[:,1]-cgdata[:,0])
cgdata[np.isnan(cgdata)] = 0.0
chrcg = {c:cgdata[cg[:,0]=='chr'+c] for c in chrom}

# scHiCluster GPU
start_time = time.time()
cluster, embedding = hicluster_gpu(network, chromsize, nc=nc)
print(time.time() - start_time)
[ARI(label, KMeans(n_clusters = nc, n_init = 200).fit(embedding[:, :ndim]).labels_) for ndim in [2,5,10,20,50]]

# scHiCluster CPU
start_time = time.time()
cluster, embedding = hicluster_cpu(network, chromsize, nc=nc, ncpus=1)
print(time.time() - start_time)
[ARI(label, KMeans(n_clusters = nc, n_init = 200).fit(embedding[:, :ndim]).labels_) for ndim in [2,5,10,20,50]]

# PCA
start_time = time.time()
cluster, embedding = raw_pca(network, chromsize, nc=nc)
print(time.time() - start_time)
[ARI(label, KMeans(n_clusters = nc, n_init = 200).fit(embedding[:, 1:(ndim+1)]).labels_) for ndim in [2,5,10,20,50]]
[ARI(label, SpectralClustering(n_clusters=nc, affinity='nearest_neighbors').fit(embedding[:, 1:(ndim+1)]).labels_) for ndim in [2,5,10,20,50]]

# Downsample reads to uniform the coverage of all the cells before PCA
start_time = time.time()
cluster, embedding = ds_pca(network, chromsize, nc=nc)
print(time.time() - start_time)
[ARI(label, KMeans(n_clusters = nc, n_init = 200).fit(embedding[:, :ndim]).labels_) for ndim in [2,5,10,20,50]]
[ARI(label, SpectralClustering(n_clusters=nc, affinity='nearest_neighbors').fit(embedding[:, :ndim]).labels_) for ndim in [2,5,10,20,50]]

# Use compartment score (PC1) of single cells
start_time = time.time()
cluster, embedding = compartment(network, chromsize, nc=nc)
print(time.time() - start_time)
[ARI(label, KMeans(n_clusters = nc, n_init = 200).fit(embedding[:, :ndim]).labels_) for ndim in [2,5,10,20,50]]
[ARI(label, SpectralClustering(n_clusters=nc, affinity='nearest_neighbors').fit(embedding[:, :ndim]).labels_) for ndim in [2,5,10,20,50]]

# Use contact-distance decay curve
start_time = time.time()
cluster, embedding = decay(network, chromsize, nc=nc)
print(time.time() - start_time)
[ARI(label, KMeans(n_clusters = nc, n_init = 200).fit(embedding[:, :ndim]).labels_) for ndim in [2,5,10,20,50]]
[ARI(label, SpectralClustering(n_clusters=nc, affinity='nearest_neighbors').fit(embedding[:, :ndim]).labels_) for ndim in [2,5,10,20,50]]

# scHiCluster without linear convolution
start_time = time.time()
cluster, embedding = hicluster_gpu(network, chromsize, nc=nc, pad=0)
print(time.time() - start_time)
[ARI(label, KMeans(n_clusters = nc, n_init = 200).fit(embedding[:, :ndim]).labels_) for ndim in [2,5,10,20,50]]
[ARI(label, SpectralClustering(n_clusters=nc, affinity='nearest_neighbors').fit(embedding[:, :ndim]).labels_) for ndim in [2,5,10,20,50]]

# scHiCluster without random walk
start_time = time.time()
cluster, embedding = hicluster_gpu(network, chromsize, nc=nc, rp=-1)
print(time.time() - start_time)
[ARI(label, KMeans(n_clusters = nc, n_init = 200).fit(embedding[:, :ndim]).labels_) for ndim in [2,5,10,20,50]]
[ARI(label, SpectralClustering(n_clusters=nc, affinity='nearest_neighbors').fit(embedding[:, :ndim]).labels_) for ndim in [2,5,10,20,50]]

# scHiCluster without keeping the top elements
start_time = time.time()
cluster, embedding = hicluster_gpu(network, chromsize, nc=nc, prct=-1)
print(time.time() - start_time)
[ARI(label, KMeans(n_clusters = nc, n_init = 200).fit(embedding[:, 1:(ndim+1)]).labels_) for ndim in [2,5,10,20,50]]
[ARI(label, SpectralClustering(n_clusters=nc, affinity='nearest_neighbors').fit(embedding[:, 1:(ndim+1)]).labels_) for ndim in [2,5,10,20,50]]
np.save('/cellar/users/zhoujt1994/projects/scHiC/'+dataset+'/embedding/1mb_pad1_rwr_real.npy', embedding)

# Random walk only
start_time = time.time()
cluster, embedding = hicluster_gpu(network, chromsize, nc=nc, pad=0, prct=-1)
print(time.time() - start_time)
[ARI(label, KMeans(n_clusters = nc, n_init = 200).fit(embedding[:, 1:(ndim+1)]).labels_) for ndim in [2,5,10,20,50]]
[ARI(label, SpectralClustering(n_clusters=nc, affinity='nearest_neighbors').fit(embedding[:, 1:(ndim+1)]).labels_) for ndim in [2,5,10,20,50]]

# Linear convolution only
start_time = time.time()
cluster, embedding = hicluster_gpu(network, chromsize, nc=nc, rp=-1, prct=-1)
print(time.time() - start_time)
[ARI(label, KMeans(n_clusters = nc, n_init = 200).fit(embedding[:, 1:(ndim+1)]).labels_) for ndim in [2,5,10,20,50]]
[ARI(label, SpectralClustering(n_clusters=nc, affinity='nearest_neighbors').fit(embedding[:, 1:(ndim+1)]).labels_) for ndim in [2,5,10,20,50]]
