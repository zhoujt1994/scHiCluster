import sys
import cooler
import numpy as np
import pandas as pd
from glob import glob
from scipy.sparse import csr_matrix, diags
from schicluster.cool.utilities import get_chrom_offsets

def load_cell_csv_to_csr(cell_path, chrom_offset, bins_df, resolution, chrom1, pos1, chrom2, pos2, min_pos_dist):
    contacts = pd.read_csv(cell_path, header=None, index_col=None, sep='\t', comment='#', 
                           usecols=[chrom1, pos1, chrom2, pos2])
    contacts = contacts[contacts[chrom1].isin(chrom_offset) & contacts[chrom2].isin(chrom_offset)]
    pos_dist = (contacts[pos1] - contacts[pos2]).abs()
    contacts = contacts[(pos_dist > min_pos_dist) |  (contacts[chrom1] != contacts[chrom2])]
    contacts['bin1_id'] = contacts[chrom1].map(chrom_offset) + (contacts[pos1] - 1) // resolution
    contacts['bin2_id'] = contacts[chrom2].map(chrom_offset) + (contacts[pos2] - 1) // resolution
    orderfilter = (contacts['bin1_id']>contacts['bin2_id'])
    contacts.loc[orderfilter, ['bin1_id', 'bin2_id']] = contacts.loc[orderfilter, ['bin2_id', 'bin1_id']].values
    
    count = contacts.groupby(['bin1_id','bin2_id'])[chrom1].count().reset_index()
    data = csr_matrix((count[chrom1].values, (count['bin1_id'].values, count['bin2_id'].values)), 
                      shape=(bins_df.shape[0], bins_df.shape[0]))
    return data

def merge_cell_raw(cell_table, chrom_size_path, output_file, resolution=5000, 
                   chrom1=1, pos1=2, chrom2=5, pos2=6, min_pos_dist=2500):
    chrom_sizes = pd.read_csv(chrom_size_path, sep='\t', index_col=0, header=None).squeeze(axis=1)
    bins_df = cooler.binnify(chrom_sizes, resolution)
    chrom_offset = get_chrom_offsets(bins_df)
    cell_list = pd.read_csv(cell_table, sep='\t', index_col=0, header=None).squeeze(axis=1)
    data = csr_matrix((bins_df.shape[0], bins_df.shape[0]))
    for xx,yy in zip(cell_list.values, cell_list.index):
        data += load_cell_csv_to_csr(cell_path=xx, chrom_offset=chrom_offset, bins_df=bins_df, resolution=resolution, 
                                     chrom1=chrom1, pos1=pos1, chrom2=chrom2, pos2=pos2, min_pos_dist=min_pos_dist)
        print(yy)

    data = data + diags(data.diagonal())
    data = data.tocoo()
    data = pd.DataFrame(np.array([data.row, data.col, data.data], dtype=int).T, 
                        columns=['bin1_id', 'bin2_id', 'count'])

    cooler.create_cooler(cool_uri=output_file, bins=bins_df, pixels=data, ordered=True)
    return
