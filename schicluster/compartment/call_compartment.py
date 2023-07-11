import subprocess
from tkinter import E
import numpy as np
import pandas as pd
from scipy.sparse import diags, csr_matrix
import cooler
import pathlib
from concurrent.futures import ProcessPoolExecutor, as_completed
import xarray as xr


def get_cpg_profile(fasta_path, hdf_output_path, cell_url=None, chrom_size_path=None, resolution=100000):
    temp_bed_path = f'{hdf_output_path}_bins.bed'
    temp_cpg_path = f'{hdf_output_path}_cpg_profile.txt'

    # count CpG per bin
    if cell_url is not None:
        cool = cooler.Cooler(cell_url)
        bins = cool.bins()[:]
    elif chrom_size_path is not None:
        chrom_sizes = pd.read_csv(chrom_size_path, sep='\t', index_col=0, header=None, squeeze=True)
        bins = cooler.binnify(chrom_sizes, resolution)
    else:
        print('ERROR : Must provide either cell_url or chrom_size_path')
    bins.to_csv(temp_bed_path, index=None, header=None, sep='\t')
    subprocess.run(
        f'bedtools nuc -fi {fasta_path} -bed {temp_bed_path} -pattern CG -C > {temp_cpg_path}',
        shell=True)

    # get CpG ratio for each bin
    result = pd.read_csv(temp_cpg_path, sep='\t')
    result['ratio'] = (result['13_user_patt_count'] /
                       (result['12_seq_len'] - result['10_num_N'])).fillna(0)
    # chrom, start, end, CpG ratio
    result = result.iloc[:, [0, 1, 2, -1]]
    result.columns = ['chrom', 'start', 'end', 'cpg_ratio']
    result.to_hdf(hdf_output_path, key='data')

    # cleanup
    subprocess.run(f'rm -f {temp_bed_path} {temp_cpg_path}', shell=True)
    return


def compartment_strength(matrix, comp, cpg_ratio):
    bin_filter = cpg_ratio > 0

    # for calculating compartment strength
    tmp = comp[bin_filter]
    a_pos = (tmp > np.percentile(tmp, 80))
    b_pos = (tmp < np.percentile(tmp, 20))
    E = matrix.tocoo()
    decay = np.array([E.diagonal(i).mean() for i in range(E.shape[0])])
    E.data = E.data / decay[np.abs(E.col - E.row)]
    E = E.tocsr()[np.ix_(bin_filter, bin_filter)]
    aa = E[np.ix_(a_pos, a_pos)].sum()
    bb = E[np.ix_(b_pos, b_pos)].sum()
    ab = E[np.ix_(a_pos, b_pos)].sum()
    scores = np.array([aa, bb, ab])
    return scores


def single_chrom_compartment(matrix, cpg_ratio, calc_strength=False):
    matrix = matrix - diags(matrix.diagonal())  # remove diag
    # matrix = matrix + matrix.T  # make full matrix
    matrix = matrix + diags((matrix.sum(axis=0).A.ravel() == 0).astype(int))
    matrix.data = matrix.data / np.repeat(matrix.sum(axis=0).A1, matrix.getnnz(axis=1))
    # weighted average of cpg density
    comp = matrix.dot(cpg_ratio.values[:, None])[:, 0]

    if calc_strength:
        scores = compartment_strength(matrix, comp, cpg_ratio)
    else:
        scores = None
    return comp, scores


def single_cell_compartment(cell_url, cpg_profile, calc_strength, output_prefix, mode,
                            resolution, chrom_sizes, chrom1, pos1, chrom2, pos2):

    if mode=='cool':
        cool = cooler.Cooler(cell_url)
        chroms = cool.chromnames
    elif mode=='tsv':
        chroms = chrom_sizes.index
        data = pd.read_csv(cell_url, sep='\t', index_col=None, header=None, comment='#')
        data = data.loc[(data[chrom1]==data[chrom2]) & data['chrom1'].isin(chroms)]
    all_comp = []
    scores = np.array([0, 0, 0])
    for chrom in chroms:
        cpg_ratio = cpg_profile.loc[cpg_profile['chrom'] == chrom, 'cpg_ratio']
        if mode=='cool':
            matrix = cool.matrix(balance=False, sparse=True).fetch(chrom)
        else:
            n_bins = (chrom_sizes.loc[chrom] // resolution) + 1
            chrfilter = (data[chrom1]==chrom)
            if chrfilter.sum()==0:
                matrix = csr_matrix((n_bins, n_bins))
            else:
                matrix = data.loc[data[chrom1]==chrfilter]
                matrix[[pos1, pos2]] = (matrix[[pos1, pos2]] - 1) // resolution
                matrix = matrix.groupby(by=[pos1, pos2])[chrom1].count().reset_index()
                matrix = csr_matrix((matrix[chrom1].astype(np.int32), (matrix[pos1], matrix[pos2])), (n_bins, n_bins))
                matrix = matrix + matrix.T
        if cpg_ratio.shape[0]!=matrix.shape[0]:
            print(f'ERROR : cpg_ratio and chrom_size have different shapes at {chrom}')
            return(0)
        comp, scores = single_chrom_compartment(matrix.tocsr(),
                                                cpg_ratio,
                                                calc_strength=calc_strength)
        if calc_strength:
            scores += scores
        all_comp.append(comp)
    all_comp = np.concatenate(all_comp)
    if calc_strength:
        np.savez(f'{output_prefix}.comp.npz', all_comp, scores)
    else:
        np.savez(f'{output_prefix}.comp.npz', all_comp)
    return


def aggregate_compartment(cell_table, temp_dir, bins, output_path, calc_strength):
    total_comp = []
    total_scores = []
    for cell_id, cell_url in cell_table.items():
        comp_path = f'{temp_dir}/{cell_id}.comp.npz'
        npz_data = np.load(comp_path)
        total_comp.append(npz_data['arr_0'])
        if calc_strength:
            total_scores.append(npz_data['arr_1'])

    total_comp = np.vstack(total_comp)
    total_comp = pd.DataFrame(total_comp,
                              index=cell_table.index,
                              columns=bins.index)
    total_comp.index.name = 'cell'
    total_comp.columns.name = 'bin'
    total_comp = xr.Dataset({'compartment': total_comp})
    total_comp.coords['bin_chrom'] = bins['chrom']
    total_comp.coords['bin_start'] = bins['start']
    total_comp.coords['bin_end'] = bins['end']

    if calc_strength:
        total_scores = np.vstack(total_scores)
        total_scores = pd.DataFrame(total_scores,
                                    index=cell_table.index,
                                    columns=['AA', 'BB', 'AB'])
        total_scores.index.name = 'cell'
        total_scores.columns.name = 'sum'
    total_comp['summary'] = xr.DataArray(total_scores)
    total_comp.to_netcdf(output_path)
    return


def multiple_cell_compartment(cell_table_path,
                              output_prefix,
                              cpg_profile_path,
                              cpu=10,
                              calc_strength=False,
                              mode='cool',
                              chrom_size_path=None,
                              resolution=100000,
                              chrom1=1,
                              pos1=2,
                              chrom2=5,
                              pos2=6):

    cell_table = pd.read_csv(cell_table_path,
                             sep='\t',
                             index_col=0,
                             header=None,
                             squeeze=True)
    print(cell_table.shape[0], 'cells to calculate.')

    temp_dir = pathlib.Path(f'{output_prefix}_compartment_temp')
    temp_dir.mkdir(exist_ok=True)
    cpg_profile = pd.read_hdf(cpg_profile_path)
    if mode=='tsv':
        if chrom_size_path is not None:
            chrom_sizes = pd.read_csv(chrom_size_path, sep='\t', index_col=0, header=None, squeeze=True)
        else:
            print('ERROR : need to provide chrom_size_path for tsv mode')
            return
    elif mode=='cool':
        chrom_sizes = None
    else:
        print('ERROR : mode need to be cool or tsv')
        return

    # calculate individual cells
    with ProcessPoolExecutor(cpu) as exe:
        future_dict = {}
        for cell_id, cell_url in cell_table.items():
            cell_prefix = f'{temp_dir}/{cell_id}'
            if pathlib.Path(f'{cell_prefix}.comp.npz').exists():
                continue
            future = exe.submit(single_cell_compartment,
                                cell_url=cell_url,
                                cpg_profile=cpg_profile,
                                output_prefix=cell_prefix,
                                calc_strength=calc_strength,
                                mode=mode,
                                chrom_sizes=chrom_sizes,
                                resolution=resolution,
                                chrom1=chrom1,
                                pos1=pos1,
                                chrom2=chrom2,
                                pos2=pos2)
            future_dict[future] = cell_id

        for future in as_completed(future_dict):
            cell_id = future_dict[future]
            print(f'{cell_id} finished.')
            future.result()

    # read bins from one of the cooler file
    # all cooler files should share the same bins
    cell_cool = cooler.Cooler(cell_table.iloc[0])
    bins = cell_cool.bins()[:]

    # aggregate boundary
    aggregate_compartment(cell_table=cell_table,
                          temp_dir=temp_dir,
                          bins=bins,
                          output_path=f'{output_prefix}.compartment.nc',
                          calc_strength=calc_strength)

    # cleanup
    subprocess.run(f'rm -rf {temp_dir}', shell=True)
    return
