import cooler
import numpy as np
import pandas as pd
from scipy.sparse import triu, csr_matrix
from concurrent.futures import ProcessPoolExecutor, as_completed

def gene_score_raw(cell_path, chrom_sizes, gene_meta, resolution, chrom1, pos1, chrom2, pos2):
    data = pd.read_csv(cell_path, sep='\t', index_col=None, header=None, comment='#')
    data = data.loc[(data[chrom1]==data[chrom2]) & data[chrom1].isin(chrom_sizes.index)]
    result = []
    for chrom in chrom_sizes.index:
        n_bins = (chrom_sizes.loc[chrom] // resolution) + 1
        chrfilter = (data[chrom1]==chrom)
        if chrfilter.sum()==0:
            D = csr_matrix((n_bins, n_bins))
        else:
            D = data.loc[data[chrom1]==chrfilter]
            D[[pos1, pos2]] = (D[[pos1, pos2]] - 1) // resolution
            D = D.groupby(by=[pos1, pos2])[chrom1].count().reset_index()
            D = csr_matrix((D[chrom1].astype(np.int32), (D[pos1], D[pos2])), (n_bins, n_bins))
        gene = gene_meta.loc[gene_meta[0]==chrom, [1,2]].values
        for xx,yy in gene:
            result.append(D[xx:(yy+1), xx:(yy+1)].sum())
    return result

def gene_score_impute(cell_path, chrom_sizes, gene_meta):
    cool = cooler.Cooler(cell_path)
    result = []
    for chrom in chrom_sizes.index:
        D = triu(cool.matrix(balance=False, sparse=True).fetch(chrom), k=1).tocsr()
        gene = gene_meta.loc[gene_meta[0]==chrom, [1,2]].values
        for xx,yy in gene:
            result.append(D[(xx-1):(yy+1), xx:(yy+2)].sum())
    return result

def gene_score(cell_table_path, gene_meta_path, resolution, output_hdf_path, chrom_size_path, 
               slop=0, cpu=10, mode='impute', chrom1=1, pos1=2, chrom2=5, pos2=6):
    chrom_sizes = pd.read_csv(chrom_size_path, sep='\t', header=None, index_col=0)
    gene_meta = pd.read_csv(gene_meta_path, sep='\t', header=None, index_col=3)
    gene_meta = gene_meta[gene_meta[0].isin(chrom_sizes.index)]
    gene_meta[1] = (gene_meta[1] - slop) // resolution
    gene_meta[2] = (gene_meta[2] + slop) // resolution
    cell_table = pd.read_csv(cell_table_path, sep='\t', header=None, index_col=None).values
    with ProcessPoolExecutor(cpu) as exe:
        future_dict = {}
        for cell, cell_path in cell_table:
            if mode=='impute':
                future = exe.submit(gene_score_impute,
                                    cell_path=cell_path,
                                    chrom_sizes=chrom_sizes,
                                    gene_meta=gene_meta)
            elif mode=='raw':
                future = exe.submit(gene_score_raw,
                                    cell_path=cell_path,
                                    chrom_sizes=chrom_sizes,
                                    gene_meta=gene_meta, 
                                    resolution=resolution, 
                                    chrom1=chrom1, 
                                    pos1=pos1, 
                                    chrom2=chrom2, 
                                    pos2=pos2)
            else:
                print("ERROR : Mode must be one of impute/raw/diag")
                return 0
            future_dict[future] = cell

        result, cell_list = [], []
        for future in as_completed(future_dict):
            cell = future_dict[future]
            print(f'{cell} finished.')
            result.append(future.result())
            cell_list.append(cell)

    result = pd.DataFrame(result, index=cell_list, columns=gene_meta.index)
    result.to_hdf(output_hdf_path, key='data', complevel=9)

