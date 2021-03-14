import time
import numpy as np
from schicluster import *
from sklearn.cluster import KMeans, SpectralClustering
from sklearn.metrics.cluster import adjusted_rand_score as ARI

mm9dim = [197195432, 181748087, 159599783, 155630120, 152537259, 149517037, 152524553, 131738871, 124076172, 129993255,
          121843856, 121257530, 120284312, 125194864, 103494974, 98319150, 95272651, 90772031, 61342430]
hg19dim = [249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747,
           135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520,
           48129895, 51304566, 155270560]

# File list and labels of dataset Ramani 2017
ctlist = ['HeLa', 'HAP1', 'GM12878', 'K562']
network = [np.loadtxt('1mb_resolution/' + ct + '/sample_list.txt', dtype=np.str) for ct in ctlist]
label = np.array([ctlist[i] for i in range(len(ctlist)) for j in range(len(network[i]))]).astype('U8')
network = np.concatenate(network)
chrom = [str(i + 1) for i in range(22)] + ['X']
chromsize = {chrom[i]: hg19dim[i] for i in range(len(chrom))}
nc = 4

# CpG content for each bin
cg = np.loadtxt('hg19/bin/hg19.1mb.bin.CpG.txt', dtype=np.str, skiprows=1, usecols=(0, 9, 11, 12))
cgdata = cg[:, 1:].astype(float)
cgdata = cgdata[:, 2] / (cgdata[:, 1] - cgdata[:, 0])
cgdata[np.isnan(cgdata)] = 0.0
chrcg = {c: cgdata[cg[:, 0] == 'chr' + c] for c in chrom}

# scHiCluster GPU
start_time = time.time()
cluster, embedding = hicluster_gpu(network, chromsize, nc=nc)
print(time.time() - start_time)
[ARI(label, KMeans(n_clusters=nc, n_init=200).fit(embedding[:, :ndim]).labels_) for ndim in [2, 5, 10, 20, 50]]

# scHiCluster CPU
start_time = time.time()
cluster, embedding = hicluster_cpu(network, chromsize, nc=nc, ncpus=5)
print(time.time() - start_time)
[ARI(label, KMeans(n_clusters=nc, n_init=200).fit(embedding[:, :ndim]).labels_) for ndim in [2, 5, 10, 20, 50]]

# PCA
start_time = time.time()
cluster, embedding = raw_pca(network, chromsize, nc=nc)
print(time.time() - start_time)
[ARI(label, KMeans(n_clusters=nc, n_init=200).fit(embedding[:, 1:(ndim + 1)]).labels_) for ndim in [2, 5, 10, 20, 50]]
[ARI(label, SpectralClustering(n_clusters=nc, affinity='nearest_neighbors').fit(embedding[:, 1:(ndim + 1)]).labels_) for
 ndim in [2, 5, 10, 20, 50]]

# Downsample reads to uniform the coverage of all the cells before PCA
start_time = time.time()
cluster, embedding = ds_pca(network, chromsize, nc=nc)
print(time.time() - start_time)
[ARI(label, KMeans(n_clusters=nc, n_init=200).fit(embedding[:, :ndim]).labels_) for ndim in [2, 5, 10, 20, 50]]
[ARI(label, SpectralClustering(n_clusters=nc, affinity='nearest_neighbors').fit(embedding[:, :ndim]).labels_) for ndim
 in [2, 5, 10, 20, 50]]

# Use compartment score (PC1) of single cells
start_time = time.time()
cluster, embedding = compartment(network, chromsize, nc=nc)
print(time.time() - start_time)
[ARI(label, KMeans(n_clusters=nc, n_init=200).fit(embedding[:, :ndim]).labels_) for ndim in [2, 5, 10, 20, 50]]
[ARI(label, SpectralClustering(n_clusters=nc, affinity='nearest_neighbors').fit(embedding[:, :ndim]).labels_) for ndim
 in [2, 5, 10, 20, 50]]

# Use contact-distance decay curve
start_time = time.time()
cluster, embedding = decay(network, chromsize, nc=nc)
print(time.time() - start_time)
[ARI(label, KMeans(n_clusters=nc, n_init=200).fit(embedding[:, :ndim]).labels_) for ndim in [2, 5, 10, 20, 50]]
[ARI(label, SpectralClustering(n_clusters=nc, affinity='nearest_neighbors').fit(embedding[:, :ndim]).labels_) for ndim
 in [2, 5, 10, 20, 50]]

# scHiCluster without linear convolution
start_time = time.time()
cluster, embedding = hicluster_gpu(network, chromsize, nc=nc, pad=0)
print(time.time() - start_time)
[ARI(label, KMeans(n_clusters=nc, n_init=200).fit(embedding[:, :ndim]).labels_) for ndim in [2, 5, 10, 20, 50]]
[ARI(label, SpectralClustering(n_clusters=nc, affinity='nearest_neighbors').fit(embedding[:, :ndim]).labels_) for ndim
 in [2, 5, 10, 20, 50]]

# scHiCluster without random walk
start_time = time.time()
cluster, embedding = hicluster_gpu(network, chromsize, nc=nc, rp=-1)
print(time.time() - start_time)
[ARI(label, KMeans(n_clusters=nc, n_init=200).fit(embedding[:, :ndim]).labels_) for ndim in [2, 5, 10, 20, 50]]
[ARI(label, SpectralClustering(n_clusters=nc, affinity='nearest_neighbors').fit(embedding[:, :ndim]).labels_) for ndim
 in [2, 5, 10, 20, 50]]

# scHiCluster without keeping the top elements
start_time = time.time()
cluster, embedding = hicluster_gpu(network, chromsize, nc=nc, prct=-1)
print(time.time() - start_time)
[ARI(label, KMeans(n_clusters=nc, n_init=200).fit(embedding[:, 1:(ndim + 1)]).labels_) for ndim in [2, 5, 10, 20, 50]]
[ARI(label, SpectralClustering(n_clusters=nc, affinity='nearest_neighbors').fit(embedding[:, 1:(ndim + 1)]).labels_) for
 ndim in [2, 5, 10, 20, 50]]
np.save('/cellar/users/zhoujt1994/projects/scHiC/' + dataset + '/embedding/1mb_pad1_rwr_real.npy', embedding)

# Random walk only
start_time = time.time()
cluster, embedding = hicluster_gpu(network, chromsize, nc=nc, pad=0, prct=-1)
print(time.time() - start_time)
[ARI(label, KMeans(n_clusters=nc, n_init=200).fit(embedding[:, 1:(ndim + 1)]).labels_) for ndim in [2, 5, 10, 20, 50]]
[ARI(label, SpectralClustering(n_clusters=nc, affinity='nearest_neighbors').fit(embedding[:, 1:(ndim + 1)]).labels_) for
 ndim in [2, 5, 10, 20, 50]]

# Linear convolution only
start_time = time.time()
cluster, embedding = hicluster_gpu(network, chromsize, nc=nc, rp=-1, prct=-1)
print(time.time() - start_time)
[ARI(label, KMeans(n_clusters=nc, n_init=200).fit(embedding[:, 1:(ndim + 1)]).labels_) for ndim in [2, 5, 10, 20, 50]]
[ARI(label, SpectralClustering(n_clusters=nc, affinity='nearest_neighbors').fit(embedding[:, 1:(ndim + 1)]).labels_) for
 ndim in [2, 5, 10, 20, 50]]
