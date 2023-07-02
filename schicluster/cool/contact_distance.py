import sys
import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed

def compute_decay(cell_name, contact_path, bins, chrom_sizes, resolution, chrom1=1, chrom2=5, pos1=2, pos2=6):
    # decay
    data = pd.read_csv(contact_path, sep='\t', header=None, index_col=None)
    data = data.loc[(data[chrom1]==data[chrom2]) & data[chrom1].isin(chrom_sizes.index)]
    hist = np.histogram(np.abs(data[pos2] - data[pos1]), bins)[0]
    # sparsity
    data[[pos1, pos2]] = data[[pos1, pos2]] // resolution
    data = data.groupby(by=[chrom1, pos1, pos2])[chrom2].count().reset_index()
    data = data.loc[data[pos1]!=data[pos2], chrom1].value_counts()
    return [pd.DataFrame(data).set_axis([cell_name], axis=1), 
            pd.DataFrame(hist, columns=[cell_name])]

def contact_distance(contact_table, chrom_size_path, resolution, output_prefix, chr1, chr2, pos1, pos2, cpu):
    chrom_sizes = pd.read_csv(chrom_size_path, sep='\t', header=None, index_col=0)
    nbins = np.floor(np.log2(chrom_sizes[1].values.max() / 2500) / 0.125)
    bins = 2500 * np.exp2(0.125 * np.arange(nbins+1))
    #dist = int(chromsize[1].min() // res + 1)
    contact_table = pd.read_csv(contact_table, sep='\t', index_col=None, header=None).values
    with ProcessPoolExecutor(cpu) as executor:
        futures = {}
        for cell_name, contact_path in contact_table.values:
            future = executor.submit(
                compute_decay,
                cell_name=cell_name,
                contact_path=contact_path,
                bins=bins,
                chrom_sizes=chrom_sizes,
                resolution=resolution,
                chrom1=chr1,
                pos1=pos1,
                chrom2=chr2,
                pos2=pos2,
            )
            futures[future] = cell_name

        sparsity, decay = {}, {}
        for future in as_completed(futures):
            cell_name = futures[future]
            xx, yy = future.result()
            sparsity[cell_name] = xx
            decay[cell_name] = yy
            print(f'{cell_name} finished')
            
    sparsity = pd.DataFrame(sparsity)[contact_table[0]].T
    decay = pd.DataFrame(decay)[contact_table[0]].T
    sparsity.to_hdf(f'{output_prefix}_chromsparsity.hdf5', key='data')
    decay.to_hdf(f'{output_prefix}_decay.hdf5', key='data')
