import pandas as pd
import numpy as np
import pathlib
import cooler
import h5py
import shutil
from concurrent.futures import ProcessPoolExecutor, as_completed

from ..cool import get_chrom_offsets
from .merge_cell_to_group import read_single_cool_chrom


def _chrom_sum_iterator(cell_urls,
                        chrom_sizes,
                        chrom_offset,
                        add_trans=False):
    """
    Iterate through the raw matrices and chromosomes of cells.

    Parameters
    ----------
    cell_urls :
        List of cell urls.
    chrom_sizes :
        Dictionary of chromosome sizes.
    chrom_offset :
        Dictionary of chromosome offsets.
    add_trans :
        If true, will also iterate all the trans combinations (different chromosomes).

    Yields
    -------
    pixel_df :
        Dataframe of pixels. Used by cooler.create_cooler to save to h5 file.
    """

    def _iter_1d(_chrom1, _chrom2):
        # sum together multiple chunks
        # first
        cell_url = cell_urls[0]
        matrix = read_single_cool_chrom(cell_url, chrom=_chrom1, chrom2=_chrom2)

        # others
        if len(cell_urls) > 1:
            for cell_url in cell_urls[1:]:
                matrix += read_single_cool_chrom(cell_url, chrom=_chrom1, chrom2=_chrom2)

        matrix = matrix.tocoo()
        _pixel_df = pd.DataFrame({
            'bin1_id': matrix.row,
            'bin2_id': matrix.col,
            'count': matrix.data
        })
        if _chrom2 is None:
            # both row and col are chrom1
            _pixel_df.iloc[:, :2] += chrom_offset[_chrom1]
        else:
            # row is chrom1, add chrom1 offset
            _pixel_df.iloc[:, 0] += chrom_offset[_chrom1]
            # col is chrom2, add chrom2 offset
            _pixel_df.iloc[:, 1] += chrom_offset[_chrom2]
        return _pixel_df

    if add_trans:
        # only iter upper triangle
        # chrom order by offset, small to large
        chroms = [k for k, v in sorted(chrom_offset.items(), key=lambda i: i[1])]
        n_chroms = len(chroms)
        for a in range(n_chroms):
            chrom1 = chroms[a]
            chrom1_dfs = []
            for b in range(a, n_chroms):
                chrom2 = chroms[b]
                pixel_df = _iter_1d(chrom1, chrom2)
                chrom1_dfs.append(pixel_df)
            chrom1_df = pd.concat(chrom1_dfs).sort_values(by=['bin1_id', 'bin2_id'])
            yield chrom1_df
    else:
        for chrom in chrom_sizes.keys():
            pixel_df = _iter_1d(chrom, None)
            yield pixel_df


def _save_single_matrix_type(cooler_path,
                             bins_df,
                             cell_urls,
                             chrom_sizes,
                             chrom_offset,
                             add_trans=False):
    """
    Save a single matrix type Cool file from merging multiple cell urls.

    Parameters
    ----------
    cooler_path :
        Path to the output cool file.
    bins_df :
        Dataframe of bins. Created from chromosome sizes and resolution.
    cell_urls :
        List of cell urls to merge.
    chrom_sizes :
        Dictionary of chromosome sizes.
    chrom_offset :
        Dictionary of chromosome offsets.
    add_trans :
        Whether add trans matrix also.
    """
    chrom_iter = _chrom_sum_iterator(cell_urls,
                                     chrom_sizes,
                                     chrom_offset,
                                     add_trans=add_trans)
    cooler.create_cooler(cool_uri=cooler_path,
                         bins=bins_df,
                         pixels=chrom_iter,
                         ordered=True,
                         dtypes={'count': np.float32})
    with h5py.File(cooler_path, 'a') as f:
        f.attrs['group_n_cells'] = len(cell_urls)
    return


def make_raw_matrix_cell_table(cell_table, resolution_str='10K'):
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
                               output_dir, add_trans=False, cpu=1):
    """
    Sum the raw matrix of cells, no normalization.

    Parameters
    ----------
    chrom_size_path :
        Path to the chrom size file. This file is used to determine chromosome names and bins.
    resolution :
        Resolution of the raw matrix.
    cell_table_path :
        Path to the cell table.
        This table should contain three columns: cell_id, cell_url, cell_group; no Header.
        The cell_id is the id of the cell in the raw matrix.
        The cell_url is the path to the raw matrix.
        The cell_group is the group of the cells.
    output_dir :
        Path to the output directory. Group cool files will be named as "<output_dir>/<cell_group>.cool".
    add_trans :
        Whether add trans matrix also.
    cpu :
        Number of CPUs to use.

    """
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
            cooler_path = output_dir / f'{cell_group}.cool'
            if cooler_path.exists():
                print(f'{cooler_path} already exists, skip.')
                continue

            cooler_temp_path = str(output_dir / f'{cell_group}.temp.cool')
            future = exe.submit(_save_single_matrix_type,
                                cooler_path=cooler_temp_path,
                                bins_df=bins_df,
                                cell_urls=cell_urls,
                                chrom_sizes=chrom_sizes,
                                chrom_offset=chrom_offset,
                                add_trans=add_trans)
            futures[future] = cell_group
        for future in as_completed(futures):
            cell_group = futures[future]
            print(f'Matrix {cell_group} generated')
            future.result()

            # move the temp cool file to the final location
            cooler_path = str(output_dir / f'{cell_group}.cool')
            cooler_temp_path = str(output_dir / f'{cell_group}.temp.cool')
            shutil.move(cooler_temp_path, cooler_path)
    return
