import sys
import cooler
import numpy as np
import pandas as pd
from glob import glob
from scipy.sparse import csr_matrix
from schicluster.cool.utilities import get_chrom_offsets

def load_cell_csv_to_csr(cell_path, chrom_offset, bins_df, resolution, chr1, pos1, chr2, pos2, min_pos_dist):
    contacts = pd.read_csv(cell_path, header=None, index_col=None, sep='\t')[[chr1, pos1, chr2, pos2]]
    contacts = contacts[contacts[chr1].isin(chrom_offset) & contacts[chr2].isin(chrom_offset)]
    pos_dist = (contacts[pos1] - contacts[pos2]).abs()
    contacts = contacts[(pos_dist > min_pos_dist) |  (contacts[chr1] != contacts[chr2])]
    contacts['bin1_id'] = contacts[chr1].map(chrom_offset) + (contacts[pos1] - 1) // resolution
    contacts['bin2_id'] = contacts[chr2].map(chrom_offset) + (contacts[pos2] - 1) // resolution
    orderfilter = (contacts['bin1_id']>contacts['bin2_id'])
    contacts.loc[orderfilter, ['bin1_id', 'bin2_id']] = contacts.loc[orderfilter, ['bin2_id', 'bin1_id']].values
    
    count = contacts.groupby(['bin1_id','bin2_id'])[chr1].count().reset_index()
    data = csr_matrix((count[chr1].values, (count['bin1_id'].values, count['bin2_id'].values)), shape=(bins_df.shape[0], bins_df.shape[0]))
    return data

def merge_cell_raw(cell_table, chrom_size_path, output_file, resolution=5000, chr1=1, pos1=2, chr2=5, pos2=6, min_pos_dist=2500):
    chrom_sizes = pd.read_csv(chrom_size_path, sep='\t', index_col=0, header=None, squeeze=True)
    bins_df = cooler.binnify(chrom_sizes, resolution)
    chrom_offset = get_chrom_offsets(bins_df)
    cell_list = pd.read_csv(cell_table, sep='\t', index_col=0)
    data = csr_matrix((bins_df.shape[0], bins_df.shape[0]))
    for xx in cell_list.index:
        data += load_cell_csv_to_csr(xx, chrom_offset, bins_df, resolution, chr1, pos1, chr2, pos2, min_pos_dist)

    data = data.tocoo()
    data = pd.DataFrame(np.array([data.row, data.col, data.data], dtype=int).T, columns=['bin1_id', 'bin2_id', 'count'])

    cooler.create_cooler(cool_uri=output_file, bins=bins_df, pixels=data, ordered=True)
    return
