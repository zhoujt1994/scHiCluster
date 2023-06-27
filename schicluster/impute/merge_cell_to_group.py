import time
import numpy as np
from scipy.sparse import csr_matrix, triu
import cooler
from ..cool.utilities import write_coo
import pandas as pd
import logging

"""
Matrix names
Q: The imputed, before normalization matrix.
"""


def read_chrom(cell_url, chrom):
    cool = cooler.Cooler(cell_url)
    matrix = triu(cool.matrix(balance=False, sparse=True).fetch(chrom))
    return matrix


def merge_cells_for_single_chromosome(cell_urls_path,
                                      chrom,
                                      output_prefix,
                                      square=False):
    cell_urls = pd.read_csv(cell_urls_path, index_col=0, header=None)[1].tolist()
    # cell_urls = cell_table['cell_url']
    n_cells = len(cell_urls)
    # get n_dims
    matrix = read_chrom(cell_urls[0], chrom)
    n_dims = matrix.shape[0]

    start_time = time.time()
    print('Merging Q (imputed, before normalization) matrix.')
    # initialize
    q_sum = csr_matrix((n_dims, n_dims), dtype=np.float32)
    q2_sum = csr_matrix((n_dims, n_dims), dtype=np.float32)
    for i, path in enumerate(cell_urls):
        data = read_chrom(path, chrom)
        q_sum += data
        if square:
            q2_sum += data.multiply(data)
    # we do not normalize by total cell numbers here, instead, normalize it in merge_group_chunks_to_group_cools
    # NO matrix_sum.data /= n_cells
    write_coo(f'{output_prefix}.Q.hdf', q_sum, chunk_size=None)
    if square:
        write_coo(f'{output_prefix}.Q2.hdf', q2_sum, chunk_size=None)
    logging.debug(f'Merge {n_cells} cells took {time.time() - start_time:.0f} seconds')
    return
