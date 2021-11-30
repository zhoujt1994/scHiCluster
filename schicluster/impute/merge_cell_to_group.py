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
                                      output_prefix):
    cell_table = pd.read_csv(cell_urls_path, index_col=0)
    cell_urls = cell_table['cell_url'].tolist()
    n_cells = len(cell_urls)
    # get n_dims
    matrix = read_chrom(cell_urls[0], chrom)
    n_dims = matrix.shape[0]

    start_time = time.time()
    print('Merging Q (imputed, before normalization) matrix.')
    # initialize
    matrix_sum = csr_matrix((n_dims, n_dims), dtype=np.float32)
    for i, path in enumerate(cell_urls):
        matrix_sum += read_chrom(path, chrom)
    # we do not normalize by total cell numbers here, instead, normalize it in merge_group_chunks_to_group_cools
    # NO matrix_sum.data /= n_cells
    write_coo(f'{output_prefix}.Q.hdf', matrix_sum, chunk_size=None)
    logging.debug(f'Merge {n_cells} cells took {time.time() - start_time:.0f} seconds')
    return
