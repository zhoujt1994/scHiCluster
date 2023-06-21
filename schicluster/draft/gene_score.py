import cooler
import numpy as np
import pandas as pd
from scipy.sparse import triu
from concurrent.futures import ProcessPoolExecutor, as_completed

def gene_score_raw(cool_path, chrom_sizes, gene):
    cool = cooler.Cooler(cool_path)
    result = []
    for chrom in chrom_sizes.index:
        D = triu(cool.matrix(balance=False, sparse=True).fetch(chrom), k=0).tocsr()
        tmp = gene.loc[gene[0]==chrom, [1,2]].values
        for xx,yy in tmp:
            result.append(D[xx:(yy+1), xx:(yy+1)].sum())
    return result

def gene_score_diag(cool_path, chrom_sizes, gene):
    cool = cooler.Cooler(cool_path)
    result = []
    for chrom in chrom_sizes.index:
        D = cool.matrix(balance=False, sparse=True).fetch(chrom).diagonal(k=0)
        tmp = gene.loc[gene[0]==chrom, [1,2]].values
        for xx,yy in tmp:
            result.append(D[xx:(yy+1)].sum())
    return result

def gene_score_impute(cool_path, chrom_sizes, gene):
    cool = cooler.Cooler(cool_path)
    result = []
    for chrom in chrom_sizes.index:
        D = triu(cool.matrix(balance=False, sparse=True).fetch(chrom), k=1).tocsr()
        tmp = gene.loc[gene[0]==chrom, [1,2]].values
        for xx,yy in tmp:
            result.append(D[(xx-1):(yy+1), xx:(yy+2)].sum())
    return result

def gene_score(cell_table, gene_meta, res, slop, output_hdf, chrom_size, cpu, mode):
    chrom_sizes = pd.read_csv(chrom_size, sep='\t', header=None, index_col=0)
    gene = pd.read_csv(gene_meta, sep='\t', header=None, index_col=3)
    gene = gene[gene[0].isin(chrom_sizes.index)]
    gene[1] = (gene[1] - slop) // res
    gene[2] = (gene[2] + slop) // res
    if mode=='impute':
        cell_table = pd.read_csv(cell_table, sep='\t', header=None, index_col=None).values
        func = gene_score_impute
    else:
        cell_table = cooler.fileops.list_coolers(cell_table)
        cell_table = np.array([[x.split('/')[-1], x] for x in cell_table])
        if mode=='raw':
            func = gene_score_raw
        elif mode=='dig':
            func = gene_score_diag
        else:
            print("ERROR : Mode must be one of impute/raw/diag")
            return(0)
    with ProcessPoolExecutor(cpu) as exe:
        future_dict = {}
        for cell, cool_path in cell_table:
            # cell_prefix = f'{temp_dir}/{cell_id}'
            # if pathlib.Path(f'{cell_prefix}.insulation.npz').exists():
            #     continue
            future = exe.submit(func,
                                cool_path,
                                chrom_sizes,
                                gene)
            future_dict[future] = cell

        result, cell_list = [], []
        for future in as_completed(future_dict):
            cell = future_dict[future]
            print(f'{cell} finished.')
            result.append(future.result())
            cell_list.append(cell)

    result = pd.DataFrame(result, index=cell_list, columns=gene.index)
    result.to_hdf(output_hdf, key='data', complevel=9)

