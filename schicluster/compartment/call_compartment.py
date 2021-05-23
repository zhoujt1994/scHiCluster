import subprocess
import numpy as np
import pandas as pd
from scipy.sparse import diags
import cooler
import pathlib
from concurrent.futures import ProcessPoolExecutor, as_completed
import xarray as xr


def get_cpg_profile(cell_url, fasta_path, hdf_output_path):
    temp_bed_path = f'{hdf_output_path}_bins.bed'
    temp_cpg_path = f'{hdf_output_path}_cpg_profile.txt'

    # count CpG per bin
    cool = cooler.Cooler(cell_url)
    bins = cool.bins()[:]
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
    matrix = matrix + matrix.T  # make full matrix
    matrix = matrix + diags((matrix.sum(axis=0).A.ravel() == 0).astype(int))

    # weighted average of cpg density
    comp = matrix.dot(cpg_ratio.values[:, None])[:, 0]

    if calc_strength:
        scores = compartment_strength(matrix, comp, cpg_ratio)
    else:
        scores = None
    return comp, scores


def single_cell_compartment(cell_url, cpg_profile_path, calc_strength,
                            output_prefix):
    cpg_profile = pd.read_hdf(cpg_profile_path)
    cool = cooler.Cooler(cell_url)

    all_comp = []
    scores = np.array([0, 0, 0])
    for chrom in cool.chromnames:
        cpg_ratio = cpg_profile.loc[cpg_profile['chrom'] == chrom, 'cpg_ratio']
        matrix = cool.matrix(balance=False, sparse=True).fetch(chrom)
        comp, scores = single_chrom_compartment(matrix,
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
                              calc_strength=False):
    cell_table = pd.read_csv(cell_table_path,
                             sep='\t',
                             index_col=0,
                             header=None,
                             squeeze=True)
    print(cell_table.shape[0], 'cells to calculate.')

    temp_dir = pathlib.Path(f'{output_prefix}_compartment_temp')
    temp_dir.mkdir(exist_ok=True)

    # calculate individual cells
    with ProcessPoolExecutor(cpu) as exe:
        future_dict = {}
        for cell_id, cell_url in cell_table.items():
            cell_prefix = f'{temp_dir}/{cell_id}'
            if pathlib.Path(f'{cell_prefix}.comp.npz').exists():
                continue
            future = exe.submit(single_cell_compartment,
                                cell_url=cell_url,
                                cpg_profile_path=cpg_profile_path,
                                output_prefix=cell_prefix,
                                calc_strength=calc_strength)
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
