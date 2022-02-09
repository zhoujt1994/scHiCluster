import sys
import numpy as np
import pandas as pd

chrom = [f'chr{i + 1}' for i in range(22)]
loop1 = pd.read_csv(sys.argv[1], sep='\t', header=None, index_col=None).sort_values(by=[0, 1, 4])
loop2 = pd.read_csv(sys.argv[2], sep='\t', header=None, index_col=None).sort_values(by=[0, 1, 4])
loop1 = loop1[loop1[0].isin(chrom)]
loop2 = loop2[loop2[0].isin(chrom)]
res = loop1.iloc[0, 2] - loop1.iloc[0, 1]
tot = 0
dist = 20000
loop0 = []
for c in chrom:
    tmp1 = loop1.loc[loop1[0] == c, [1, 4]].values
    tmp2 = loop2.loc[loop2[0] == c, [1, 4]].values
    j, tmp = 0, 0
    for i in range(len(tmp1)):
        while j <= (len(tmp2) - 1) and tmp2[j, 0] < (tmp1[i, 0] - dist):
            j += 1
        k = j
        flag = 1
        while k <= (len(tmp2) - 1) and tmp2[k, 0] <= (tmp1[i, 0] + dist):
            if (tmp1[i, 1] - dist) <= tmp2[k, 1] <= (tmp1[i, 1] + dist):
                tmp += 1
                flag = 0
                break
            k += 1
        if flag:
            loop0.append([c, tmp1[i, 0], tmp1[i, 0] + res, c, tmp1[i, 1], tmp1[i, 1] + res])
    tot += tmp

print(tot, loop1.shape[0], tot / loop1.shape[0])
if len(sys.argv) == 3:
    np.savetxt('.'.join(sys.argv[1].split('.')[:-1]) + '.specific.bedpe',
               loop0, fmt='%s', delimiter='\t')
