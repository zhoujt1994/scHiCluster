import numpy as np
import pandas as pd
import cooler
from sklearn.decomposition import TruncatedSVD
import joblib
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
import pathlib


def make_idx(n_dim, dist, resolution):
    idx = np.triu_indices(n_dim, k=1)
    idx_filter = np.array([(yy - xx) < (dist / resolution + 1)
                           for xx, yy in zip(idx[0], idx[1])])
    idx = (idx[0][idx_filter], idx[1][idx_filter])
    return idx


def make_chrom_matrix(cell_table,
                      chrom,
                      nbins,
                      output_path,
                      scale_factor,
                      dist,
                      resolution):
    idx = make_idx(nbins, dist, resolution)
    shape = (cell_table.size, idx[0].size)
    # read data
    chrom_matrix = np.zeros(shape=shape, dtype='float32')
    for i, (_, cell_url) in enumerate(cell_table.items()):
        cool = cooler.Cooler(cell_url)
        matrix = cool.matrix(balance=False, sparse=False).fetch(chrom)
        # each row of chrom_matrix is a 1D cell-chrom matrix
        chrom_matrix[i, :] = matrix[idx].ravel()
    chrom_matrix *= scale_factor
    np.savez(output_path, chrom_matrix)
    return


def svd(input_path, dim, output_prefix, save_model=True, norm_sig=True):
    chrom_matrix = np.load(input_path)['arr_0']
    dim = min(dim, chrom_matrix.shape[0] - 1, chrom_matrix.shape[1] - 1)
    model = TruncatedSVD(n_components=dim, algorithm='arpack')
    decomp = model.fit_transform(chrom_matrix)

    if norm_sig:
        decomp = decomp[:, model.singular_values_ > 0]
        singular_values = model.singular_values_[model.singular_values_ > 0]
        decomp /= singular_values[None, :]

    np.savez(f'{output_prefix}_decomp.npz', decomp)
    if save_model:
        joblib.dump(model, f'{output_prefix}_SVD.lib')
    return


def embedding(cell_table_path,
              output_dir,
              dim=50,
              dist=1000000,
              resolution=100000,
              scale_factor=100000,
              norm_sig=True,
              cpu=1,
              save_model=False,
              save_raw=True):
    cell_table = pd.read_csv(cell_table_path,
                             sep='\t',
                             index_col=0,
                             header=None,
                             squeeze=True)

    output_dir = pathlib.Path(output_dir).absolute()
    raw_dir = output_dir / 'raw'
    raw_dir.mkdir(exist_ok=True)

    first_cool = cooler.Cooler(cell_table.iloc[0])
    chroms = first_cool.chromnames
    chrom_bin_counts = first_cool.bins()[:]['chrom'].value_counts()
    # remove small chroms
    chroms = [chrom for chrom in chroms if chrom_bin_counts[chrom] > 2]

    # prepare raw chromosome 1D matrix
    with ProcessPoolExecutor(cpu) as exe:
        futures = {}
        for chrom in chroms:
            nbins = chrom_bin_counts[chrom]
            output_path = raw_dir / f'{chrom}.npz'
            future = exe.submit(make_chrom_matrix,
                                cell_table=cell_table,
                                chrom=chrom,
                                nbins=nbins,
                                output_path=output_path,
                                scale_factor=scale_factor,
                                dist=dist,
                                resolution=resolution)
            futures[future] = chrom

        for future in as_completed(futures):
            chrom = futures[future]
            print(f'{chrom} matrix generated')
            future.result()

    # SVD on each chromosome
    decomp_dir = output_dir / 'decomp'
    decomp_dir.mkdir(exist_ok=True)
    with ProcessPoolExecutor(cpu) as exe:
        futures = {}
        for chrom in chroms:
            chrom_raw_path = raw_dir / f'{chrom}.npz'
            output_prefix = decomp_dir / chrom
            future = exe.submit(svd,
                                input_path=chrom_raw_path,
                                dim=dim,
                                output_prefix=output_prefix,
                                save_model=save_model,
                                norm_sig=norm_sig)
            futures[future] = chrom

        for future in as_completed(futures):
            chrom = futures[future]
            print(f'{chrom} SVD generated')
            future.result()

    # concatenate all chromosome decomp matrix and do another SVD
    total_data = []
    for chrom in chroms:
        decomp_path = decomp_dir / f'{chrom}_decomp.npz'
        data = np.load(str(decomp_path))['arr_0']
        total_data.append(data)
    total_data = np.concatenate(total_data, axis=1)
    total_data_path = decomp_dir / f'total_chrom_decomp_concat.npz'
    np.savez(total_data_path, total_data)

    # final SVD
    output_prefix = decomp_dir / f'total'
    svd(input_path=total_data_path,
        dim=dim,
        output_prefix=output_prefix,
        save_model=save_model,
        norm_sig=norm_sig)

    # clean single chrom
    for chrom in chroms:
        decomp_path = decomp_dir / f'{chrom}_decomp.npz'
        subprocess.run(['rm', '-f', str(decomp_path)])
    if not save_raw:
        subprocess.run(['rm', '-rf', str(raw_dir)])
    return
