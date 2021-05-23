import pathlib
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
import pandas as pd
from rpy2.rinterface_lib.embedded import RRuntimeError
from rpy2.robjects import r, pandas2ri, numpy2ri
from rpy2.robjects.packages import importr, isinstalled
from rpy2.robjects.vectors import StrVector
import cooler
from schicluster.cool.utilities import get_chrom_offsets
from scipy.sparse import csr_matrix, save_npz, load_npz, vstack
import anndata
import xarray as xr
import subprocess
import schicluster

pandas2ri.activate()
numpy2ri.activate()

PACKAGE_DIR = pathlib.Path(schicluster.__path__[0])


# check if some package is installed
def install_r_package(name):
    if not isinstalled(name):
        utils = importr('utils')
        utils.chooseCRANmirror(ind=1)
        utils.install_packages(StrVector([name]))


def domain_df_to_boundary(cool, total_results, resolution):
    # aggr to 1D boundary record,
    # 1 means either start or end of a domain,
    # 2 means both start and end of two domains
    bins = cool.bins()[:]
    chrom_offset = get_chrom_offsets(bins)
    use_results = total_results[total_results['name'] == 'domain']
    bin_loc = use_results.iloc[:, [1, 2]] // resolution
    bin_loc += use_results['chrom'].map(chrom_offset).values[:, None]
    bin_loc_data = np.zeros(bins.shape[0], dtype=np.int16)
    bin_loc_data[bin_loc['chromStart']] += 1
    bin_loc_data[bin_loc['chromEnd']] += 1
    bin_loc_data = csr_matrix(bin_loc_data)
    return bin_loc_data


def single_chrom_calculate_insulation_score(matrix, window_size=10):
    w = window_size
    score = np.ones(matrix.shape[0])
    for i in range(1, matrix.shape[0]):
        if i < w:
            intra = (matrix[:i, :i].sum() +
                     matrix[i:(i + w), i:(i + w)].sum()) / (i * (i + 1) / 2 +
                                                            (w * (w + 1) / 2))
            inter = matrix[:i, i:(i + w)].sum() / (i * (i + w))
        else:
            intra = (matrix[(i - w):i, (i - w):i].sum() +
                     matrix[i:(i + w), i:(i + w)].sum()) / (w * (w + 1))
            inter = matrix[(i - w):i, i:(i + w)].sum() / (w * w)
        score[i] = inter / (inter + intra)
    return score


def call_domain_and_insulation(cell_url,
                               output_prefix,
                               resolution=25000,
                               window_size=10):
    r.source(str(PACKAGE_DIR / 'domain/TopDom.R'))
    pandas2ri.activate()
    numpy2ri.activate()

    def run_top_dom(matrix, bins):
        j = matrix.indices + 1
        p = matrix.indptr
        x = matrix.data
        result = r.RunTopDom(j, p, x, bins, window_size)
        result = pd.DataFrame(result)
        return result

    cool = cooler.Cooler(cell_url)
    total_domain_results = []
    total_insulation_score = []
    for chrom in cool.chromnames:
        matrix = cool.matrix(balance=False, sparse=True).fetch(chrom).tocsc()
        bins = cool.bins().fetch(chrom).reset_index(drop=True)
        bins.columns = ["chr", "from.coord", "to.coord"]
        if (matrix.nnz < (matrix.shape[0] * matrix.shape[1] * 0.001)) or (matrix.shape[0] < 10):
            # skip if matrix is too less info or too small
            pass
        else:
            try:
                total_domain_results.append(run_top_dom(matrix, bins))
            except RRuntimeError:
                print('Got R error at', cell_url, chrom, matrix.shape, matrix.data.size, bins.shape)
        total_insulation_score.append(
            single_chrom_calculate_insulation_score(matrix, window_size))
    if len(total_domain_results) == 0:
        total_domain_results = pd.DataFrame([], columns=['chrom', 'chromStart', 'chromEnd', 'name'])
    else:
        total_domain_results = pd.concat(total_domain_results).reset_index(
            drop=True)
    domain_boundary = domain_df_to_boundary(cool, total_domain_results,
                                            resolution)
    insulation_score = np.concatenate(total_insulation_score)

    # save
    save_npz(f'{output_prefix}.boundary.npz', domain_boundary)
    np.savez(f'{output_prefix}.insulation.npz', insulation_score)
    return


def aggregate_boundary(cell_table, temp_dir, bins, output_path):
    total_boundary = []

    for cell_id, cell_url in cell_table.items():
        boundary_path = f'{temp_dir}/{cell_id}.boundary.npz'
        total_boundary.append(load_npz(boundary_path))
    total_boundary = vstack(total_boundary)
    adata = anndata.AnnData(X=total_boundary,
                            obs=pd.DataFrame([], index=cell_table.index.tolist()),
                            var=bins)
    adata.write_h5ad(output_path)
    return


def aggregate_insulation(cell_table, temp_dir, bins, output_path):
    total_insulation = []
    for cell_id, cell_url in cell_table.items():
        insulation_path = f'{temp_dir}/{cell_id}.insulation.npz'
        total_insulation.append(np.load(insulation_path)['arr_0'])
    total_insulation = np.vstack(total_insulation)
    total_insulation = pd.DataFrame(total_insulation,
                                    index=cell_table.index,
                                    columns=bins.index)
    total_insulation.index.name = 'cell'
    total_insulation.columns.name = 'bin'
    total_insulation = xr.DataArray(total_insulation)
    total_insulation.coords['bin_chrom'] = bins['chrom']
    total_insulation.coords['bin_start'] = bins['start']
    total_insulation.coords['bin_end'] = bins['end']
    total_insulation.to_netcdf(output_path)
    return


def multiple_call_domain_and_insulation(cell_table_path,
                                        output_prefix,
                                        resolution=25000,
                                        window_size=10,
                                        cpu=10):
    # install R package Matrix
    install_r_package('Matrix')

    cell_table = pd.read_csv(cell_table_path,
                             sep='\t',
                             index_col=0,
                             header=None,
                             squeeze=True)
    print(cell_table.shape[0], 'cells to calculate.')

    temp_dir = pathlib.Path(f'{output_prefix}_domain_temp')
    temp_dir.mkdir(exist_ok=True)

    # calculate individual cells
    with ProcessPoolExecutor(cpu) as exe:
        future_dict = {}
        for cell_id, cell_url in cell_table.items():
            cell_prefix = f'{temp_dir}/{cell_id}'
            if pathlib.Path(f'{cell_prefix}.insulation.npz').exists():
                continue
            future = exe.submit(call_domain_and_insulation,
                                cell_url,
                                cell_prefix,
                                resolution=resolution,
                                window_size=window_size)
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
    aggregate_boundary(cell_table=cell_table,
                       temp_dir=temp_dir,
                       bins=bins,
                       output_path=f'{output_prefix}.boundary.h5ad')

    # aggregate insulation
    aggregate_insulation(cell_table=cell_table,
                         temp_dir=temp_dir,
                         bins=bins,
                         output_path=f'{output_prefix}.insulation.nc')

    # cleanup
    subprocess.run(f'rm -rf {temp_dir}', shell=True)
    return
