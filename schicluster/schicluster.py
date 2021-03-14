import os
import sys
import time
import numpy as np
import torch
import torch.nn.functional as F
from multiprocessing import Pool
from scipy.sparse import csr_matrix
from scipy.stats import chi2_contingency
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans


def neighbor_ave_gpu(A, pad):
    if pad == 0:
        return torch.from_numpy(A).float().cuda()
    ll = pad * 2 + 1
    conv_filter = torch.ones(1, 1, ll, ll).cuda()
    B = F.conv2d(torch.from_numpy(A[None, None, :, :]).float().cuda(), conv_filter, padding=pad * 2)
    return (B[0, 0, pad:-pad, pad:-pad] / float(ll * ll))


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
    cell, c, ngene, pad, rp = args
    D = np.loadtxt(cell + '_chr' + c + '.txt')
    A = csr_matrix((D[:, 2], (D[:, 0], D[:, 1])), shape=(ngene, ngene)).toarray()
    A = np.log2(A + A.T + 1)
    A = neighbor_ave_gpu(A, pad)
    if rp == -1:
        Q = A[:]
    else:
        Q = random_walk_gpu(A, rp)
    return Q.reshape(ngene * ngene)


def hicluster_gpu(network, chromsize, nc, res=1000000, pad=1, rp=0.5, prct=20, ndim=20):
    matrix = []
    for i, c in enumerate(chromsize):
        ngene = int(chromsize[c] / res) + 1
        start_time = time.time()
        Q_concat = torch.zeros(len(network), ngene * ngene).float().cuda()
        for j, cell in enumerate(network):
            Q_concat[j] = impute_gpu([cell, c, ngene, pad, rp])
        Q_concat = Q_concat.cpu().numpy()
        if prct > -1:
            thres = np.percentile(Q_concat, 100 - prct, axis=1)
        Q_concat = (Q_concat > thres[:, None])
        end_time = time.time()
        print('Load and impute chromosome', c, 'take', end_time - start_time, 'seconds')
        ndim = int(min(Q_concat.shape) * 0.2) - 1
        # U, S, V = torch.svd(Q_concat, some=True)
        # R_reduce = torch.mm(U[:, :ndim], torch.diag(S[:ndim])).cuda().numpy()
        pca = PCA(n_components=ndim)
        R_reduce = pca.fit_transform(Q_concat)
        matrix.append(R_reduce)
        print(c)
    matrix = np.concatenate(matrix, axis=1)
    pca = PCA(n_components=min(matrix.shape) - 1)
    matrix_reduce = pca.fit_transform(matrix)
    kmeans = KMeans(n_clusters=nc, n_init=200).fit(matrix_reduce[:, :ndim])
    return kmeans.labels_, matrix_reduce


def neighbor_ave_cpu(A, pad):
    if pad == 0:
        return A
    ngene, _ = A.shape
    ll = pad * 2 + 1
    B, C, D, E = [np.zeros((ngene + ll, ngene + ll)) for i in range(4)]
    B[(pad + 1):(pad + ngene + 1), (pad + 1):(pad + ngene + 1)] = A[:]
    F = B.cumsum(axis=0).cumsum(axis=1)
    C[ll:, ll:] = F[:-ll, :-ll]
    D[ll:, :] = F[:-ll, :]
    E[:, ll:] = F[:, :-ll]
    return (np.around(F + C - D - E, decimals=8)[ll:, ll:] / float(ll * ll))


def random_walk_cpu(A, rp):
    ngene, _ = A.shape
    A = A - np.diag(np.diag(A))
    A = A + np.diag(np.sum(A, axis=0) == 0)
    P = np.divide(A, np.sum(A, axis=0))
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
    cell, c, ngene, pad, rp = args
    D = np.loadtxt(cell + '_chr' + c + '.txt')
    A = csr_matrix((D[:, 2], (D[:, 0], D[:, 1])), shape=(ngene, ngene)).toarray()
    A = np.log2(A + A.T + 1)
    A = neighbor_ave_cpu(A, pad)
    if rp == -1:
        Q = A[:]
    else:
        Q = random_walk_cpu(A, rp)
    return [cell, Q.reshape(ngene * ngene)]


def hicluster_cpu(network, chromsize, nc, res=1000000, pad=1, rp=0.5, prct=20, ndim=20, ncpus=10):
    matrix = []
    for i, c in enumerate(chromsize):
        ngene = int(chromsize[c] / res) + 1
        start_time = time.time()
        paras = [[cell, c, ngene, pad, rp] for cell in network]
        p = Pool(ncpus)
        result = p.map(impute_cpu, paras)
        p.close()
        index = {x[0]: j for j, x in enumerate(result)}
        Q_concat = np.array([result[index[x]][1] for x in network])
        if prct > -1:
            thres = np.percentile(Q_concat, 100 - prct, axis=1)
            Q_concat = (Q_concat > thres[:, None])
        end_time = time.time()
        print('Load and impute chromosome', c, 'take', end_time - start_time, 'seconds')
        ndim = int(min(Q_concat.shape) * 0.2) - 1
        pca = PCA(n_components=ndim)
        R_reduce = pca.fit_transform(Q_concat)
        matrix.append(R_reduce)
        print(c)
    matrix = np.concatenate(matrix, axis=1)
    pca = PCA(n_components=min(matrix.shape) - 1)
    matrix_reduce = pca.fit_transform(matrix)
    kmeans = KMeans(n_clusters=nc, n_init=200).fit(matrix_reduce[:, :ndim])
    return kmeans.labels_, matrix_reduce


def raw_pca(network, chromsize, nc, res=1000000, ndim=20):
    matrix = []
    for i, c in enumerate(chromsize):
        ngene = int(chromsize[c] / res) + 1
        start_time = time.time()
        uptri = np.triu_indices(ngene, 1)
        A_concat = np.zeros((len(label), len(uptri[0]))).astype(float)
        j = 0
        for cell in network:
            D = np.loadtxt(cell + '_chr' + c + '.txt')
            A = csr_matrix((D[:, 2], (D[:, 0], D[:, 1])), shape=(ngene, ngene)).toarray()
            A = np.log2(A + A.T + 1)
            A_concat[j, :] = A[uptri]
            j += 1
        end_time = time.time()
        print('Load and impute chromosome', c, 'take', end_time - start_time, 'seconds')
        pca = PCA(n_components=min(A_concat.shape) - 1)
        A_reduce = pca.fit_transform(A_concat)
        matrix.append(A_reduce)
        print(c)
    matrix = np.concatenate(matrix, axis=1)
    pca = PCA(n_components=min(matrix.shape) - 1)
    matrix_reduce = pca.fit_transform(matrix)
    kmeans = KMeans(n_clusters=nc, n_init=200).fit(matrix_reduce[:, 1:(ndim + 1)])
    return kmeans.labels_, matrix_reduce


def ds_pca(network, chromsize, nc, res=1000000, ndim=20):
    matrix = []
    for i, c in enumerate(chromsize):
        ngene = int(chromsize[c] / res) + 1
        start_time = time.time()
        uptri = np.triu_indices(ngene, 1)
        A_concat = np.zeros((len(label), len(uptri[0]))).astype(float)
        j = 0
        tot = int(np.min([np.sum(np.loadtxt(cell + '_chr' + c + '.txt')[:, 2]) for cell in network]))
        for cell in network:
            D = np.loadtxt(cell + '_chr' + c + '.txt')
            A = csr_matrix((ngene, ngene)).toarray().astype(float)
            idx = np.concatenate([[j for i in range(int(D[j, 2]))] for j in range(len(D))])
            shu = np.arange(len(idx))
            np.random.shuffle(shu)
            for x in idx[shu[:tot]]:
                A[int(D[x, 0]), int(D[x, 1])] += 1
            A = np.log2(A + A.T + 1)
            A_concat[j, :] = A[uptri]
            j += 1
        end_time = time.time()
        print('Load and impute chromosome', c, 'take', end_time - start_time, 'seconds')
        pca = PCA(n_components=min(A_concat.shape) - 1)
        A_reduce = pca.fit_transform(A_concat)
        matrix.append(A_reduce)
        print(c)
    matrix = np.concatenate(matrix, axis=1)
    pca = PCA(n_components=min(matrix.shape) - 1)
    matrix_reduce = pca.fit_transform(matrix)
    kmeans = KMeans(n_clusters=nc, n_init=200).fit(matrix_reduce[:, :ndim])
    return kmeans.labels_, matrix_reduce


def comp_comp(O, cg):
    ngene, _ = Q.shape
    N = torch.zeros(ngene, ngene).float().cuda()
    for k in range(ngene):
        diagsum = torch.sum(torch.diag(O, k)) / (ngene - k)
        N = N + torch.diag(torch.diag(O, k) / diagsum, k)
    D = torch.diag(torch.diag(N))
    N = N + torch.t(N) - D
    N[N != N] = 0.0
    C = corrcoef(N).cpu().numpy()
    pca = PCA(n_components=2)
    pc = pca.fit_transform(C)
    corr = np.abs([np.corrcoef(cg, pc[:, i])[0, 1] for i in range(2)])
    # if corr[0] > corr[1]:
    # 	k = 0
    # else:
    # 	k = 1
    k = 0
    pc, corr = pc[:, k], corr[k]
    comp = torch.mean(N[pc > 0][:, pc < 0]).item() * 2 / (
            torch.mean(N[pc > 0][:, pc > 0]).item() + torch.mean(N[pc < 0][:, pc < 0]).item())
    if np.corrcoef(cg, pc)[0, 1] < 0:
        return [C, -pc, corr, comp]
    else:
        return [C, pc, corr, comp]


def compartment(network, chromsize, nc, res=1000000, ndim=20):
    global chrcg
    pca = PCA(n_components=2)
    matrix = []
    for i, c in enumerate(chrom):
        ngene = int(chromsize[c] / res) + 1
        start_time = time.time()
        comp = np.zeros((len(label), ngene)).astype(float)
        j = 0
        for cell in network:
            D = np.loadtxt(cell + '_chr' + c + '.txt')
            A = csr_matrix((D[:, 2], (D[:, 0], D[:, 1])), shape=(ngene, ngene)).toarray()
            A = np.log2(A + A.T + 1)
            A = neighbor_ave_gpu(A, 0)
            comp[j] = comp_comp(A, chrcg[c])
            j += 1
        end_time = time.time()
        print('Load and random walk take', end_time - start_time, 'seconds')
        matrix.append(comp)
        print(c)
    matrix = np.concatenate(matrix, axis=1)
    pca = PCA(n_components=min(matrix.shape) - 1)
    matrix_reduce = pca.fit_transform(matrix)
    kmeans = KMeans(n_clusters=nc, n_init=200).fit(matrix_reduce[:, :ndim])
    return kmeans.labels_, matrix_reduce


def decay(network, chromsize, nc, res=1000000, ndim=20):
    matrix = []
    for i, c in enumerate(chromsize):
        ngene = int(chromsize[c] / res) + 1
        start_time = time.time()
        dec = np.zeros((len(label), ngene - 1)).astype(float)
        j = 0
        for j, cell in enumerate(network):
            D = np.loadtxt(cell + '_chr' + c + '.txt')
            A = csr_matrix((D[:, 2], (D[:, 0], D[:, 1])), shape=(ngene, ngene)).toarray()
            tmp = np.array([np.sum(np.diag(A, k)) for k in range(1, ngene)])
            dec[j] = tmp / np.sum(tmp)
        end_time = time.time()
        print('Load and random walk take', end_time - start_time, 'seconds')
        matrix.append(dec)
        print(c)
    matrix = np.concatenate(matrix, axis=1)
    pca = PCA(n_components=min(matrix.shape) - 1)
    matrix_reduce = pca.fit_transform(matrix)
    kmeans = KMeans(n_clusters=nc, n_init=200).fit(matrix_reduce[:, :ndim])
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
    np.savetxt(cell + '_chr' + c + '.topdommatrix', C, fmt='%s', delimiter='\t')
    return


def output_sparse(cell, c, Q, res):
    idx = np.where(Q > 0)
    C = np.concatenate(((idx[0] * res)[:, None], (idx[1] * res)[:, None], Q[idx][:, None]), axis=1)
    np.savetxt(cell + '_chr' + c + '.sparsematrix', C, fmt='%d\t%d\t%s', delimiter='\n')
    return


def diff_dom(args):
    domlist, cluster, ctlist, res, c, chromsize = args
    start_time = time.time()
    ctdict = {x: i for i, x in enumerate(ctlist)}
    ngene = int(chromsize[c] / res) + 1
    bound = np.zeros((len(cluster), ngene))
    for i, cell in enumerate(domlist):
        domain = np.loadtxt(cell, dtype=np.str)
        C = domain[domain[:, -1] == 'domain', :2].astype(int) / res
        if len(C) > int(chromsize[c] / 10000000) + 1:
            bound[i, np.array(list(set(C.flatten())), dtype=int)] += 1
    cellfilter = (np.sum(bound, axis=1) > 0)
    count = np.array([np.sum(np.logical_and(cluster == k, cellfilter)) for k in ctlist])
    prob = np.array([np.sum(bound[cluster == k], axis=0) for k in ctlist])
    ptmp, btmp = [], []
    for i in range(ngene):
        contig = [[prob[j, i], count[j] - prob[j, i]] for j in range(ctlist)]
        if np.sum(np.sum(contig, axis=0) == 0) == 0:
            ptmp.append(chi2_contingency(contig)[1])
            btmp.append('\t'.join([c, str(i * res), str(i * res + res)]))
    print(c, time.time() - start_time)
    return [bound, prob / count[:, None], btmp, ptmp]


def chrflankpv(args):
    bins, prob, cluster, ctlist, res, c, chromsize, w = args
    start_time = time.time()
    sel = (np.sum(prob, axis=1) > 0)
    tmpclust = cluster[sel]
    prob = prob[sel]
    flankpv = []
    for x in bins:
        x = x.split('\t')
        mid = int(x[1]) // res
        tmpcount = np.max(prob[:, (mid - w):(mid + w + 1)], axis=1)
        contig = np.array(
            [[np.sum(tmpcount[tmpclust == j]), np.sum(tmpclust == j) - np.sum(tmpcount[tmpclust == j])] for j in
             ctlist])
        r = contig[:, 0] / np.sum(contig, axis=1)
        if np.sum(contig == 0) == 0:
            flankpv.append(x + [chi2_contingency(contig)[1]] + r)
    flankpv = np.array(flankpv)
    print(c, time.time() - start_time)
    return [flankpv[:, 4:].astype(float), flankpv[:, 3].astype(float), flankpv[:, :3]]


def filter_bins(prob, fdr, fdr_cutoff=0.05, diff_cutoff=0, max_cutoff=0, min_cutoff=1):
    maxp = np.max(prob, axis=0)
    minp = np.min(prob, axis=0)
    ans = np.logical_and(maxp > max_cutoff, minp < min_cutoff)
    ans = np.logical_and(ans, maxp - minp > diff_cutoff)
    ans = np.logical_and(ans, fdr < fdr_cutoff)
    return ans
