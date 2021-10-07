import logging
import pathlib

import cooler
import numpy as np
import pandas as pd


def get_chrom_offsets(bins_df):
    chrom_offset = {chrom: bins_df[bins_df['chrom'] == chrom].index[0]
                    for chrom in bins_df['chrom'].cat.categories}
    return chrom_offset


def write_coo(path, matrix, chunk_size=5000000):
    """
    Write chromosome contacts as chunked pixel df (cooler input)
    """
    matrix = matrix.tocoo(copy=False)
    df = pd.DataFrame({'bin1_id': matrix.row, 'bin2_id': matrix.col, 'count': matrix.data})
    with pd.HDFStore(path, complib='zlib', complevel=3) as hdf:
        if chunk_size is None:
            # no chunk
            hdf['c0'] = df
        else:
            for i, chunk_start in enumerate(range(0, df.shape[0], chunk_size)):
                hdf[f'c{i}'] = df.iloc[chunk_start:chunk_start + chunk_size]
    return


def chrom_iterator(input_dir, chrom_order, chrom_offset, chrom_wildcard='{chrom}.hdf'):
    for chrom in chrom_order:
        chrom_file = chrom_wildcard.format(chrom=chrom)
        output_path = f'{input_dir}/{chrom_file}'
        if not pathlib.Path(output_path).exists():
            continue
        with pd.HDFStore(output_path, 'r') as hdf:
            logging.debug(chrom)
            keys = {int(i[2:]): i for i in hdf.keys()}
            for i in sorted(keys.keys()):
                key = keys[i]
                chunk = hdf[key]
                chunk.iloc[:, :2] += chrom_offset[chrom]
                yield chunk


def aggregate_chromosomes(chrom_size_path,
                          resolution,
                          input_dir,
                          output_path,
                          chrom_wildcard='{chrom}.hdf'):
    chrom_sizes = cooler.read_chromsizes(chrom_size_path, all_names = True)
    bins_df = cooler.binnify(chrom_sizes, resolution)
    chrom_offset = get_chrom_offsets(bins_df)

    cooler.create_cooler(cool_uri=output_path,
                         bins=bins_df,
                         pixels=chrom_iterator(input_dir=input_dir,
                                               chrom_order=bins_df['chrom'].unique(),
                                               chrom_offset=chrom_offset,
                                               chrom_wildcard=chrom_wildcard),
                         ordered=True,
                         dtypes={'count': np.float32})
    return


def cell_chunk(cell_url, chrom_sizes, chunk=50000000):
    cell_cool = cooler.Cooler(cell_url)
    chunk_df = cooler.binnify(chrom_sizes, chunk)
    for _, row in chunk_df.iterrows():
        chrom, start, end = row
        region = f'{chrom}:{start}-{end}'
        # this fetch selected a rectangle region, row is region, col is whole chrom
        data = cell_cool.matrix(balance=False, as_pixels=True).fetch(region, chrom)
        yield data


def aggregate_cells(output_path, cell_dir, chrom_size_path, resolution):
    chrom_sizes = cooler.read_chromsizes(chrom_size_path, all_names=True)
    bins_df = cooler.binnify(chrom_sizes, resolution)

    cell_pixel_dict = {
        cell_url.name.split('.')[0]: cell_chunk(str(cell_url), chrom_sizes=chrom_sizes)
        for cell_url in pathlib.Path(cell_dir).glob('*cool')
    }
    cooler.create_scool(output_path,
                        bins=bins_df,
                        cell_name_pixels_dict=cell_pixel_dict,
                        ordered=True,
                        mode='a')
    return
