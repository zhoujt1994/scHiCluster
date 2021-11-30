import pandas as pd
import numpy as np
import pathlib
import cooler
import h5py
from concurrent.futures import ProcessPoolExecutor, as_completed

from ..cool import get_chrom_offsets
from .merge_cell_to_group import read_single_cool_chrom


def chrom_sum_iterator(cell_urls,
                       chrom_sizes,
                       chrom_offset):
    for chrom in chrom_sizes.keys():
        # sum together multiple chunks
        # first
        cell_url = cell_urls[0]
        matrix = read_single_cool_chrom(cell_url, chrom)
        # others
        for cell_url in cell_urls[1:]:
            matrix += read_single_cool_chrom(cell_url, chrom)

        matrix = matrix.tocoo()
        pixel_df = pd.DataFrame({
            'bin1_id': matrix.row,
            'bin2_id': matrix.col,
            'count': matrix.data
        })
        pixel_df.iloc[:, :2] += chrom_offset[chrom]
        yield pixel_df


def save_single_matrix_type(cooler_path,
                            bins_df,
                            cell_urls,
                            chrom_sizes,
                            chrom_offset):
    chrom_iter = chrom_sum_iterator(cell_urls,
                                    chrom_sizes,
                                    chrom_offset)
    cooler.create_cooler(cool_uri=cooler_path,
                         bins=bins_df,
                         pixels=chrom_iter,
                         ordered=True,
                         dtypes={'count': np.float32})
    with h5py.File(cooler_path, 'a') as f:
        f.attrs['group_n_cells'] = len(cell_urls)
    return


def make_raw_matrix_cell_table(cell_table_path, resolution_str='10K'):
    cell_table = pd.read_csv(cell_table_path,
                             header=None,
                             sep='\t',
                             index_col=0,
                             names=['cell_id', 'cell_url', 'cell_group'])

    # get all the raw matrix cell url automatically if the dir path is the default
    scool_dirs = cell_table['cell_url'].apply(lambda i: '/'.join(i.split('/')[:-4])).unique()
    cell_urls = {}
    scool_file_pattern = f'raw/*.{resolution_str}.scool'
    for scool_dir in scool_dirs:
        for scool_path in pathlib.Path(scool_dir).glob(scool_file_pattern):
            with h5py.File(scool_path, 'r') as _cool:
                cell_ids = list(_cool['cells'].keys())
                for cell_id in cell_ids:
                    cell_urls[cell_id] = f'{scool_path}::/cells/{cell_id}'
    raw_url_series = pd.Series(cell_urls)

    # delete old url
    del cell_table['cell_url']
    # add new url
    cell_table['cell_url'] = raw_url_series

    na_cells = cell_table['cell_url'].isna().sum()
    if na_cells > 0:
        raise ValueError(f'{na_cells} cells do not have raw matrix.')

    cell_table = cell_table[['cell_url', 'cell_group']]
    return cell_table


def merge_raw_scool_by_cluster(chrom_size_path, resolution, cell_table_path,
                               output_dir, cpu=1):
    """Sum the raw matrix of cells, no normalization."""
    # determine chunk dirs for the group:
    output_dir = pathlib.Path(output_dir).absolute()
    output_dir.mkdir(exist_ok=True)
    cell_table = pd.read_csv(cell_table_path,
                             sep='\t',
                             index_col=0,
                             header=None,
                             names=['cell_id', 'cell_url', 'cell_group'])

    chrom_sizes = cooler.read_chromsizes(chrom_size_path, all_names=True)
    bins_df = cooler.binnify(chrom_sizes, resolution)
    chrom_offset = get_chrom_offsets(bins_df)

    with ProcessPoolExecutor(cpu) as exe:
        futures = {}
        for cell_group, sub_df in cell_table.groupby('cell_group'):
            cell_urls = sub_df['cell_url'].tolist()
            cooler_path = str(output_dir / f'{cell_group}.cool')
            future = exe.submit(save_single_matrix_type,
                                cooler_path=cooler_path,
                                bins_df=bins_df,
                                cell_urls=cell_urls,
                                chrom_sizes=chrom_sizes,
                                chrom_offset=chrom_offset)
            futures[future] = cell_group
        for future in as_completed(futures):
            cell_group = futures[future]
            print(f'Matrix {cell_group} generated')
            future.result()
    return
