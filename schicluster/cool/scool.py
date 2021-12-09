import re
import subprocess
import warnings
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed

import cooler
import pandas.errors

from .utilities import get_chrom_offsets
import pandas as pd
from cooler import create_scool
from .remove_blacklist import filter_contacts


def generate_scool_batch_data(cell_path_dict,
                              resolution,
                              chrom_offset,
                              chrom_size_path,
                              blacklist_1d_path,
                              blacklist_2d_path,
                              remove_duplicates,
                              blacklist_resolution,
                              output_path,
                              chr1=1,
                              chr2=5,
                              pos1=2,
                              pos2=6,
                              min_pos_dist=2500):
    def single_cell_pixel(_path):
        try:
            contacts = filter_contacts(_path,
                                       chrom_size_path=chrom_size_path,
                                       blacklist_1d_path=blacklist_1d_path,
                                       blacklist_2d_path=blacklist_2d_path,
                                       remove_duplicates=remove_duplicates,
                                       resolution_2d=blacklist_resolution,
                                       chrom1=chr1,
                                       pos1=pos1,
                                       chrom2=chr2,
                                       pos2=pos2)
        except pandas.errors.EmptyDataError:
            # empty contacts file
            return pd.DataFrame([], columns=['bin1_id', 'bin2_id', 'count'])
        pos_dist = (contacts[pos1] - contacts[pos2]).abs()

        # filter
        contacts = contacts[contacts[chr1].isin(chrom_offset)
                            & contacts[chr2].isin(chrom_offset)
                            & (contacts[pos1] > 0)
                            & (contacts[pos2] > 0)
                            & (pos_dist > min_pos_dist)]

        # calculate pixel
        bin1_id = contacts[chr1].map(chrom_offset) + (contacts[pos1] - 1) // resolution
        bin2_id = contacts[chr2].map(chrom_offset) + (contacts[pos2] - 1) // resolution
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
                                     blacklist_1d_path,
                                     blacklist_2d_path,
                                     remove_duplicates,
                                     blacklist_resolution,
                                     chr1=1,
                                     chr2=5,
                                     pos1=2,
                                     pos2=6,
                                     min_pos_dist=2500,
                                     batch_n=20,
                                     cpu=1):
    # parse chromosome sizes, prepare bin_df
    chrom_sizes = pd.read_csv(
        chrom_size_path,
        sep='\t',
        index_col=0,
        header=None, squeeze=True)
    bins_df = cooler.binnify(chrom_sizes, resolution)
    chrom_offset = get_chrom_offsets(bins_df)

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
                           output_path=batch_output,
                           chr1=chr1, chr2=chr2,
                           pos1=pos1, pos2=pos2,
                           min_pos_dist=min_pos_dist,
                           chrom_size_path=chrom_size_path,
                           blacklist_1d_path=blacklist_1d_path,
                           blacklist_2d_path=blacklist_2d_path,
                           remove_duplicates=remove_duplicates,
                           blacklist_resolution=blacklist_resolution)
            futures[f] = batch_output

        for future in as_completed(futures):
            # batch finished
            batch_output = futures[future]
            future.result()

            # dump batch result into cool
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
                   blacklist_1d_path=None,
                   blacklist_2d_path=None,
                   blacklist_resolution=10000,
                   remove_duplicates=True,
                   chr1=1,
                   chr2=5,
                   pos1=2,
                   pos2=6,
                   min_pos_dist=2500,
                   cpu=1,
                   batch_n=50):
    """
    Generate single-resolution cool files from single-cell contact files recorded in contacts_table

    Parameters
    ----------
    contacts_table
        tab-separated table containing tow columns, 1) cell id, 2) cell contact file path (juicer-pre format)
        No header
    output_prefix
        Output prefix of the cool files. Output path will be {output_prefix}.{resolution_str}.cool
    chrom_size_path
        Path to the chromosome size file, this file should be in UCSC chrom.sizes format. We will use this file as
        the reference to create matrix. It is recommended to remove small contigs or chrM from this file.
    resolutions
        Resolutions to generate the matrix. Each resolution will be stored in a separate file.
    blacklist_1d_path
        Path to blacklist region BED file, such as ENCODE blacklist.
        Either side of the contact overlapping with a blacklist region will be removed.
    blacklist_2d_path
        Path to blacklist region pair BEDPE file.
        Both side of the contact overlapping with the same blacklist region pair will be removed.
    blacklist_resolution
        Resolution in bps when consider the 2D blacklist region pairs.
    remove_duplicates
        If true, will remove duplicated contacts based on [chr1, pos1, chr2, pos2] values
    chr1
        0 based index of chr1 column.
    chr2
        0 based index of chr2 column.
    pos1
        0 based index of pos1 column.
    pos2
        0 based index of pos2 column.
    min_pos_dist
        Minimum distance for a fragment to be considered.
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
                                         cpu=cpu,
                                         chr1=chr1, chr2=chr2,
                                         pos1=pos1, pos2=pos2,
                                         min_pos_dist=min_pos_dist,
                                         blacklist_1d_path=blacklist_1d_path,
                                         blacklist_2d_path=blacklist_2d_path,
                                         remove_duplicates=remove_duplicates,
                                         blacklist_resolution=blacklist_resolution)
        print('Finished', output_path)
    return
