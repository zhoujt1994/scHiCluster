import pandas as pd
from cooler import create_scool
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
import warnings
import subprocess
import re


def get_chrom_offsets(chrom_sizes_series, resolution):
    cur_offset = 0
    chrom_offset = {}
    for chrom, length in chrom_sizes_series.items():
        chrom_offset[chrom] = cur_offset
        n_bins = length // resolution + 1 * (length % resolution != 0)
        cur_offset += n_bins
    return chrom_offset


def generate_bins_df_from_chrom_sizes(chrom_sizes, resolution):
    bins_dfs = []
    for chrom, length in chrom_sizes.items():
        chrom_bins_df = pd.DataFrame({
            'start':
                list(range(0, length, resolution)),
            'end':
                list(range(resolution, length + resolution, resolution))
        })
        chrom_bins_df.iloc[-1, 1] = length
        chrom_bins_df['chrom'] = chrom
        bins_dfs.append(chrom_bins_df)
    bins_df = pd.concat(bins_dfs)
    bins_df = bins_df[['chrom', 'start',
                       'end']].sort_values(['chrom',
                                            'start']).reset_index(drop=True)
    return bins_df


def generate_scool_batch_data(cell_path_dict, resolution, chrom_offset,
                              output_path):
    def single_cell_pixel(path):
        contacts = pd.read_csv(path, sep='\t', header=None)

        # filter
        contacts = contacts[contacts[1].isin(chrom_offset)
                            & contacts[5].isin(chrom_offset)
                            & (contacts[2] > 0)
                            & (contacts[6] > 0)]

        # calculate pixel
        bin1_id = contacts[1].map(chrom_offset) + (contacts[2] -
                                                   1) // resolution
        bin2_id = contacts[5].map(chrom_offset) + (contacts[6] -
                                                   1) // resolution
        counter = Counter()
        for x, y in zip(bin1_id, bin2_id):
            if x > y:
                x, y = y, x
            counter[(x, y)] += 1
        pixel_df = pd.Series(counter).reset_index()
        pixel_df.columns = pd.Index(['bin1_id', 'bin2_id', 'count'])
        pixel_df = pixel_df.sort_values(['bin1_id',
                                         'bin2_id']).reset_index(drop=True)
        return pixel_df

    with warnings.catch_warnings():
        # ignore the
        warnings.simplefilter("ignore")
        with pd.HDFStore(output_path) as hdf:
            for cell_id, path in cell_path_dict.items():
                hdf[cell_id] = single_cell_pixel(path)
    return


def generate_scool_single_resolution(cell_path_dict,
                                     chrom_size_path,
                                     resolution,
                                     output_path,
                                     batch_n=20,
                                     cpu=1):
    # parse chromosome sizes, prepare bin_df
    chrom_sizes = pd.read_csv(chrom_size_path,
                              header=None,
                              index_col=0,
                              sep='\t',
                              squeeze=True).sort_index()
    chrom_offset = get_chrom_offsets(chrom_sizes, resolution)
    bins_df = generate_bins_df_from_chrom_sizes(chrom_sizes, resolution)

    chunk_dicts = defaultdict(dict)
    for i, (cell, path) in enumerate(cell_path_dict.items()):
        chunk_dicts[i // batch_n][cell] = path

    with ProcessPoolExecutor(cpu) as exe:
        futures = {}
        for batch, cell_path_dict in chunk_dicts.items():
            batch_output = output_path + f'_{batch}'
            f = exe.submit(generate_scool_batch_data,
                           cell_path_dict=cell_path_dict,
                           resolution=resolution,
                           chrom_offset=chrom_offset,
                           output_path=batch_output)
            futures[f] = batch_output

        for future in as_completed(futures):
            # batch finished
            batch_output = futures[future]
            future.result()

            # dump batch result into scool
            cell_pixel_dict = {}
            with pd.HDFStore(batch_output, mode='r') as hdf:
                for cell_id in hdf.keys():
                    cell_id = cell_id[1:]  # remove '/'
                    cell_pixel_dict[cell_id] = hdf[cell_id]
            create_scool(output_path,
                         bins=bins_df,
                         cell_name_pixels_dict=cell_pixel_dict,
                         ordered=True,
                         mode='a')
            subprocess.run(['rm', '-f', batch_output], check=True)
    return


def generate_scool(contacts_table,
                   output_prefix,
                   chrom_size_path,
                   resolutions,
                   cpu=1,
                   batch_n=50):
    """
    Generate single-resolution scool files from single-cell contact files recorded in contacts_table

    Parameters
    ----------
    contacts_table
        tab-separated table containing tow columns, 1) cell id, 2) cell contact file path (juicer-pre format)
        No header
    output_prefix
        Output prefix of the scool files. Output path will be {output_prefix}.{resolution_str}.scool
    chrom_size_path
        Path to the chromosome size file, this file should be in UCSC chrom.sizes format. We will use this file as
        the reference to create matrix. It is recommended to remove small contigs or chrM from this file.
    resolutions
        Resolutions to generate the matrix. Each resolution will be stored in a separate file.
    cpu
        number of cpus to parallel.
    batch_n
        number of cells to deal with in each cpu process.

    Returns
    -------

    """
    # read cell paths
    cell_path_dict = pd.read_csv(contacts_table,
                                 header=None,
                                 index_col=0,
                                 squeeze=True,
                                 sep='\t').to_dict()
    print(f'{len(cell_path_dict)} cells to process.')

    for resolution in resolutions:
        resolution_str = str(resolution)
        resolution_str = re.sub(r'000000$', 'M', resolution_str)
        resolution_str = re.sub(r'000$', 'K', resolution_str)
        output_path = f'{output_prefix}.{resolution_str}.scool'

        print('Generating', output_path)
        generate_scool_single_resolution(cell_path_dict=cell_path_dict,
                                         chrom_size_path=chrom_size_path,
                                         resolution=resolution,
                                         output_path=output_path,
                                         batch_n=batch_n,
                                         cpu=cpu)
        print('Finished', output_path)
    return
