import time
import h5py
import cooler
import pathlib
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, load_npz, triu
from ..cool import write_coo, get_chrom_offsets
from concurrent.futures import ProcessPoolExecutor, as_completed

"""
Matrix names

"""


def merge_cells_for_single_chromosome(output_dir,
                                      output_prefix,
                                      merge_type='E'):
    """
    Merge cell's E and T matrix to group matrices, sum only, not normalized by n_cells yet
    E: Matrix normalized by global diagonal backgrounds, calculated from loop_bkg
    E2: Sum of square of E, used to calculate global t statistics
    T: Matrix normalized by global diagonal and local backgrounds, then minus E (T is the delta matrix),
    calculated from loop_bkg
    T2: Sum of square of T, used to calculate local t statistics
    """
    start_time = time.time()
    # get cell paths
    cell_paths = [str(p) for p in pathlib.Path(output_dir).glob(f'*.{merge_type}.npz')]
    n_cells = len(cell_paths)
    # get n_dims
    matrix = load_npz(cell_paths[0])
    n_dims = matrix.shape[0]
    # initialize
    e_sum = csr_matrix((n_dims, n_dims), dtype=np.float32)
    e2_sum = csr_matrix((n_dims, n_dims), dtype=np.float32)
    for i, chrom_path in enumerate(cell_paths):
        matrix = load_npz(chrom_path)
        e_sum += matrix
        e2_sum += matrix.multiply(matrix)
    write_coo(f'{output_prefix}.{merge_type}.hdf', e_sum, chunk_size=None)
    write_coo(f'{output_prefix}.{merge_type}2.hdf', e2_sum, chunk_size=None)
    print(f'Merge {n_cells} cells took {time.time() - start_time:.0f} seconds')
    return


def read_single_cool_chrom(cool_path, chrom, chrom2=None):
    # Used in chrom_sum_iterator, return the sum according to group_n_cells
    # Also used in merge_raw_matrix and merge_group
    # Output chrom matrix
    cool = cooler.Cooler(str(cool_path))
    selector = cool.matrix(balance=False, sparse=True)
    if chrom2 is None:
        matrix = triu(selector.fetch(chrom))
    else:
        if chrom == chrom2:
            matrix = triu(selector.fetch(chrom, chrom2))
        else:
            matrix = selector.fetch(chrom, chrom2)
    with h5py.File(cool_path, 'r') as f:
        if 'group_n_cells' in f.attrs:
            matrix.data *= f.attrs['group_n_cells']
    return matrix


def chrom_sum_iterator(input_cool_list,
                       chrom_sizes,
                       chrom_offset, 
                       total_cells):
    # Used in save_single_matrix_type, return the average over cells
    # total_cells need to be provided
    # Output chrom df
    for chrom in chrom_sizes.keys():
        cool_path = input_cool_list[0]
        matrix = read_single_cool_chrom(cool_path, chrom)
        for cool_path in input_cool_list[1:]:
            matrix += read_single_cool_chrom(cool_path, chrom)
        matrix = matrix.tocoo()
        pixel_df = pd.DataFrame({
            'bin1_id': matrix.row,
            'bin2_id': matrix.col,
            'count': matrix.data
        })
        pixel_df.iloc[:, :2] += chrom_offset[chrom]
        pixel_df.iloc[:, -1] /= total_cells
        yield pixel_df


def save_single_matrix_type(input_cool_list, 
                            output_cool,
                            bins_df,
                            chrom_sizes,
                            chrom_offset,
                            total_cells):
    # Used by merge_group_chunks_to_group_cools and merge_cool
    # total_cells need to be provided
    # Output cool
    chrom_iter = chrom_sum_iterator(input_cool_list,
                                    chrom_sizes,
                                    chrom_offset,
                                    total_cells)
    cooler.create_cooler(cool_uri=output_cool,
                         bins=bins_df,
                         pixels=chrom_iter,
                         ordered=True,
                         dtypes={'count': np.float32})
    with h5py.File(output_cool, 'a') as f:
        f.attrs['group_n_cells'] = total_cells
    return
 

def merge_cool(input_cool_tsv_file, output_cool):
    # Input could be cool files of single cell or average over cells
    # Output is average over cells
    # total_cell is counted over cools according to group_n_cells, otherwise 1
    input_cool_list = pd.read_csv(input_cool_tsv_file, header=None, index_col=None).squeeze(axis=1).tolist()
    input_cool_list = [str(pathlib.Path(cool).absolute()) for cool in input_cool_list]
    
    cool = cooler.Cooler(input_cool_list[0])
    bins_df = cool.bins()[["chrom", "start", "end"]][:]
    chrom_sizes = cool.chromsizes[:]
    chrom_offset = get_chrom_offsets(bins_df)
    total_cells = 0
    for cool_path in input_cool_list:
        with h5py.File(cool_path, 'r') as f:
            if 'group_n_cells' in f.attrs:
                total_cells += f.attrs['group_n_cells']
            else:
                total_cells += 1

    save_single_matrix_type(input_cool_list, 
                            output_cool,
                            bins_df,
                            chrom_sizes,
                            chrom_offset,
                            total_cells)
    return


def merge_group_chunks_to_group_cools(chrom_size_path,
                                      resolution,
                                      group,
                                      output_dir,
                                      matrix_types=('E', 'E2', 'T', 'T2', 'Q', 'Q2')):
    # Input is sum over cells per chunk
    # Output is average over cells of all chunks
    # total_cell is counted over cell_table.csv per chunk

    # determine chunk dirs for the group:
    output_dir = pathlib.Path(output_dir).absolute()
    group_dir = output_dir / group
    group_dir.mkdir(exist_ok=True)

    chunk_dirs = list(output_dir.glob(f'{group}_chunk*'))

    # count total cells
    total_cells = 0
    for chunk_dir in chunk_dirs:
        total_cells += pd.read_csv(chunk_dir / 'cell_table.csv', index_col=0, header=None).shape[0]

    chrom_sizes = cooler.read_chromsizes(chrom_size_path, all_names=True)
    bins_df = cooler.binnify(chrom_sizes, resolution)
    chrom_offset = get_chrom_offsets(bins_df)

    with ProcessPoolExecutor(6) as exe:
        futures = {}
        for matrix_type in matrix_types:
            output_cool = str(group_dir / f'{group}.{matrix_type}.cool')
            input_cool_list = [chunk_dir / f'{matrix_type}.cool' for chunk_dir in chunk_dirs]
            future = exe.submit(save_single_matrix_type,
                                input_cool_list=input_cool_list,
                                output_cool=output_cool,
                                bins_df=bins_df,
                                chrom_sizes=chrom_sizes,
                                chrom_offset=chrom_offset,
                                total_cells=total_cells)
            futures[future] = matrix_type
        for future in as_completed(futures):
            matrix_type = futures[future]
            print(f'Matrix {matrix_type} generated')
            future.result()
    return


def merge_group_to_bigger_group_cools(chrom_size_path,
                                      resolution,
                                      group,
                                      output_dir,
                                      group_list,
                                      shuffle,
                                      matrix_types=('E', 'E2', 'T', 'T2', 'Q', 'Q2')):
    """
    Sum all the group average cool files,
    and finally divide the total number of cells to
    get a group cell number normalized cool file in the end.
    """
    # determine chunk dirs for the group:
    output_dir = pathlib.Path(output_dir).absolute()
    group_dir = output_dir / group
    group_dir.mkdir(exist_ok=True, parents=True)

    group_list = [pathlib.Path(xx) for xx in group_list]
    
    # count total cells
    total_cells = 0
    for chunk_dir in group_list:
        total_cells += pd.read_csv(chunk_dir / 'cell_table.csv', index_col=0, header=None).shape[0]
    
    if shuffle:
        group_list = [xx / 'shuffle/' for xx in group_list]

    chrom_sizes = cooler.read_chromsizes(chrom_size_path, all_names=True)
    bins_df = cooler.binnify(chrom_sizes, resolution)
    chrom_offset = get_chrom_offsets(bins_df)

    with ProcessPoolExecutor(6) as exe:
        futures = {}
        for matrix_type in matrix_types:
            output_cool = str(group_dir / f'{group}.{matrix_type}.cool')
            input_cool_list = [list(group_dir.glob(f'*/*.{matrix_type}.cool'))[0] for group_dir in group_list]
            future = exe.submit(save_single_matrix_type,
                                input_cool_list=input_cool_list,
                                output_cool=output_cool,
                                bins_df=bins_df,
                                chrom_sizes=chrom_sizes,
                                chrom_offset=chrom_offset,
                                total_cells=total_cells)
            futures[future] = matrix_type
        for future in as_completed(futures):
            matrix_type = futures[future]
            print(f'Matrix {matrix_type} generated')
            future.result()
    return
