import pathlib
import time
import numpy as np
from scipy.sparse import csr_matrix, load_npz, triu
from ..cool import write_coo, get_chrom_offsets
from .merge_cell_to_group import read_single_cool_chrom
import cooler
import pandas as pd
import h5py
from concurrent.futures import ProcessPoolExecutor, as_completed

def chrom_ave_iterator(chunk_dirs,
                       chrom_sizes,
                       chrom_offset,
                       matrix_type,
                       total_cells):
    print(f'Reading matrix {matrix_type}')
    for chrom in chrom_sizes.keys():
        # sum together multiple chunks
        # first
        cool_path = list(chunk_dirs[0].glob(f'*/*.{matrix_type}.cool'))[0]
        with h5py.File(cool_path, 'a') as f:
            n_cells = f.attrs['group_n_cells']
        matrix = read_single_cool_chrom(cool_path, chrom) * n_cells
        # others
        for chunk_dir in chunk_dirs[1:]:
            cool_path = list(chunk_dir.glob(f'*/*.{matrix_type}.cool'))[0]
            with h5py.File(cool_path, 'a') as f:
                n_cells = f.attrs['group_n_cells']
            matrix += read_single_cool_chrom(cool_path, chrom) * n_cells

        matrix = matrix.tocoo()
        pixel_df = pd.DataFrame({
            'bin1_id': matrix.row,
            'bin2_id': matrix.col,
            'count': matrix.data
        })
        pixel_df.iloc[:, :2] += chrom_offset[chrom]
        # All the chunk
        pixel_df.iloc[:, -1] /= total_cells
        yield pixel_df


def save_single_matrix_type(cooler_path,
                            bins_df,
                            chunk_dirs,
                            chrom_sizes,
                            chrom_offset, 
                            matrix_type,
                            total_cells):
    chrom_iter = chrom_ave_iterator(chunk_dirs,
                                    chrom_sizes,
                                    chrom_offset,
                                    matrix_type,
                                    total_cells)
    cooler.create_cooler(cool_uri=cooler_path,
                         bins=bins_df,
                         pixels=chrom_iter,
                         ordered=True,
                         dtypes={'count': np.float32})
    with h5py.File(cooler_path, 'a') as f:
        f.attrs['group_n_cells'] = total_cells
    return


def merge_group_to_bigger_group_cools(chrom_size_path,
                                      resolution,
                                      group,
                                      output_dir,
                                      group_list,
                                      shuffle,
                                      matrix_types=('E', 'E2', 'T', 'T2', 'Q')):
    """
    Sum all the chunk sum cool files,
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
        total_cells += pd.read_csv(chunk_dir / 'cell_table.tsv', index_col=0).shape[0]
    
    if shuffle:
        group_list = [xx / 'shuffle/' for xx in group_list]

    chrom_sizes = cooler.read_chromsizes(chrom_size_path, all_names=True)
    bins_df = cooler.binnify(chrom_sizes, resolution)
    chrom_offset = get_chrom_offsets(bins_df)

    with ProcessPoolExecutor(5) as exe:
        futures = {}
        for matrix_type in matrix_types:
            cooler_path = str(group_dir / f'{group}.{matrix_type}.cool')
            future = exe.submit(save_single_matrix_type,
                                cooler_path=cooler_path,
                                bins_df=bins_df,
                                chunk_dirs=group_list,
                                chrom_sizes=chrom_sizes,
                                chrom_offset=chrom_offset,
                                matrix_type=matrix_type,
                                total_cells=total_cells)
            futures[future] = matrix_type
        for future in as_completed(futures):
            matrix_type = futures[future]
            print(f'Matrix {matrix_type} generated')
            future.result()
    return


