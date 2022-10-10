import pathlib
import shutil
from concurrent.futures import ProcessPoolExecutor, as_completed

import cooler
import numpy as np
import pandas as pd
import xarray as xr
import zarr
from cooler.util import read_chromsizes, binnify
from numcodecs import Blosc

SMALL_SAMPLE_CHUNK = 1
COMPRESSOR_C_LEVEL = 3


class CoolDSSingleMatrixWriter:
    def __init__(self,
                 path,
                 cool_table_path,
                 value_types,
                 chrom_sizes_path,
                 chrom1,
                 chrom2=None,
                 mode='w',
                 cooler_bin_size=10000,
                 bin_chunk_size=510,
                 sample_chunk_size=50,
                 data_dtype='float32',
                 cpu=1):
        """
        Write a single chrom1-by-chrom2 matrix to CoolDS zarr.

        Parameters
        ----------
        path :
            Path to the zarr dir.
        cool_table_path :
            Path to the cool table with four columns: sample, value_type, path, cool_type
        value_types :
            Dict of cool types and their value types.
        chrom_sizes_path :
            Path to the chrom sizes file.
        chrom1 :
            Chrom1 name.
        chrom2 :
            Chrom2 name. If None, chrom1 will be used.
        mode :
            Mode to open the zarr.
        cooler_bin_size :
            Cooler bin size.
        bin_chunk_size :
            Chunk size of the bin1 and bin2 dimensions.
        sample_chunk_size :
            Chunk size of the sample dimension.
        data_dtype :
            Data type of the matrix.
        cpu :
            Number of CPUs to use.
        """
        self.path = path
        self.root = zarr.open(path, mode=mode)
        self.log_dir_path = None
        self.mode = mode
        self.value_types = value_types
        self.cool_tables, self.sample_ids = self._read_cool_table(cool_table_path)
        self.n_sample = self.sample_ids.size
        self.n_cpu = cpu

        self.bin_chunk_size = bin_chunk_size
        self.sample_chunk_size = sample_chunk_size
        self.data_dtype = data_dtype

        self.cooler_bin_size = cooler_bin_size
        self.chrom_sizes, self.chrom_bins, self.n_bins = self._read_chrom_info(chrom_sizes_path)

        self.chrom1 = chrom1
        if chrom2 is None:
            self.chrom2 = chrom1
        else:
            self.chrom2 = chrom2

        self.chrom1_size = self.chrom_sizes[self.chrom1]
        self.chrom2_size = self.chrom_sizes[self.chrom2]
        self.chrom1_bins = self.chrom_bins[self.chrom_bins['chrom'] == self.chrom1]
        self.chrom2_bins = self.chrom_bins[self.chrom_bins['chrom'] == self.chrom2]
        self.chrom1_n_bins = self.chrom1_bins.shape[0]
        self.chrom2_n_bins = self.chrom2_bins.shape[0]

        self.execute()

    def _read_cool_table(self, cool_table_path):
        value_types = self.value_types
        # get cool tables for each cool type
        cool_paths = pd.read_csv(cool_table_path,
                                 header=None,
                                 names=['sample', 'value_type', 'path', 'cool_type'])

        _cool_tables = {}
        sample_ids = []
        for cool_type, sub_paths in cool_paths.groupby('cool_type'):
            cool_table = sub_paths.pivot(index='sample',
                                         columns='value_type',
                                         values='path')

            assert cool_type in value_types, f'cool_type "{cool_type}" not in value_types provided by user.'
            value_type_list = value_types[cool_type]
            assert all([value in cool_table.columns for value in value_type_list
                        ]), f'Some {cool_type} values in {value_type_list} missing in cool_table'
            assert cool_table.isna().values.sum(
            ) == 0, f'{cool_type} cool path table has some missing paths.'

            _cool_tables[cool_type] = cool_table
            sample_ids += cool_table.index.tolist()

        sample_ids = pd.Index(set(sample_ids))
        cool_tables = {}
        for cool_type, table in _cool_tables.items():
            assert table.index.size == sample_ids.size, 'sample id not consistent between cool types'
            cool_tables[cool_type] = table.loc[sample_ids].copy()

        return cool_tables, sample_ids

    def _read_chrom_info(self, chrom_sizes_path):
        # chrom length and bins info
        chrom_sizes = read_chromsizes(chrom_sizes_path)
        chrom_bins = binnify(chrom_sizes, self.cooler_bin_size)
        n_bins = chrom_bins.shape[0]
        return chrom_sizes, chrom_bins, n_bins

    def _add_root_attrs(self):
        self.root.attrs['chrom1'] = self.chrom1
        self.root.attrs['chrom2'] = self.chrom2
        self.root.attrs['chrom1_size'] = self.chrom1_size
        self.root.attrs['chrom2_size'] = self.chrom2_size
        self.root.attrs['chrom1_n_bins'] = self.chrom1_n_bins
        self.root.attrs['chrom2_n_bins'] = self.chrom2_n_bins
        self.root.attrs['cooler_bin_size'] = self.cooler_bin_size
        self.root.attrs['bin_chunk_size'] = self.bin_chunk_size
        self.root.attrs['sample_chunk_size'] = self.sample_chunk_size

    def _init_zarr(self):
        root = zarr.open(self.path, mode=self.mode)

        # add log dir
        self.log_dir_path = pathlib.Path(self.path) / '.log'
        self.log_dir_path.mkdir(exist_ok=True, parents=True)

        # prepare flag
        success_flag = self.log_dir_path / "_PREPARE_SUCCESS"
        if success_flag.exists():
            print("Zarr structure already initiated. Skipping...")

        # create sample_id
        sample_id = root.require_dataset("sample_id",
                                         shape=(self.n_sample,),
                                         chunks=(self.n_sample,),
                                         dtype="<U256")
        sample_id[:] = list(self.sample_ids)
        sample_id.attrs["_ARRAY_DIMENSIONS"] = "sample_id"

        # create empty da
        for cool_type, value_type_list in self.value_types.items():
            # value type index
            n_value_type = len(value_type_list)
            value_dim = f"{cool_type}_value_type"
            value_type_da = root.require_dataset(value_dim,
                                                 shape=(n_value_type,),
                                                 chunks=(n_value_type,),
                                                 dtype="<U10")
            value_type_da[:] = list(value_type_list)
            value_type_da.attrs["_ARRAY_DIMENSIONS"] = [value_dim]

            # data
            z = root.require_dataset(
                cool_type,
                shape=(self.chrom1_n_bins, self.chrom2_n_bins, self.n_sample, n_value_type),
                chunks=(self.bin_chunk_size, self.bin_chunk_size, self.sample_chunk_size, 1),
                dtype=self.data_dtype
            )
            z.attrs["_ARRAY_DIMENSIONS"] = ["bin1", "bin2", "sample_id", value_dim]

        # add root attrs
        self._add_root_attrs()

        # consolidate metadata
        zarr.consolidate_metadata(self.path)

        # touch success flag
        success_flag.touch()
        return

    def _save_cool_to_temp_zarr(self, cool_paths, output_path, cool_type, value_type, triu):
        da_list = []
        for sample, cool_path in cool_paths.items():
            cool = cooler.Cooler(cool_path)
            matrix = cool.matrix(balance=False, sparse=False).fetch(self.chrom1, self.chrom2).astype(self.data_dtype)
            assert matrix.shape == (self.chrom1_n_bins, self.chrom2_n_bins)

            # keep only the upper triangle of the matrix
            if triu:
                matrix = np.triu(matrix)

            # da.dims {'bin1': n_bins, 'bin2': n_bins, 'sample_id': 1, 'value_type': 1}
            da = xr.DataArray(matrix, dims=['bin1', 'bin2']).chunk({
                'bin1': self.bin_chunk_size,
                'bin2': self.bin_chunk_size
            }).expand_dims(dim={
                'sample_id': [sample],
                f'{cool_type}_value_type': [value_type]
            }).transpose('bin1', 'bin2', 'sample_id', f'{cool_type}_value_type')
            da_list.append(da)
        da = xr.concat(da_list, dim='sample_id')
        ds = xr.Dataset({cool_type: da})
        ds.to_zarr(output_path, mode='w')
        return

    def _save_small_sample_chunks(self):
        with ProcessPoolExecutor(self.n_cpu) as executor:
            futures = {}
            for cool_type, cool_table in self.cool_tables.items():
                for value_type, sample_cool_path in cool_table.items():
                    for chunk_start in range(0, self.n_sample, SMALL_SAMPLE_CHUNK):
                        sample_chunk = self.sample_ids[chunk_start:chunk_start + SMALL_SAMPLE_CHUNK]
                        cool_paths = sample_cool_path.loc[sample_chunk]
                        output_path = f'{self.chrom1}_{self.chrom2}_{cool_type}_{value_type}_{chunk_start}_temp.zarr'
                        if self.chrom2 is None or self.chrom1 == self.chrom2:
                            # save upper triangle only if matrix is cis contacts
                            triu = True
                        else:
                            triu = False
                        future = executor.submit(self._save_cool_to_temp_zarr,
                                                 cool_paths=cool_paths,
                                                 output_path=output_path,
                                                 cool_type=cool_type,
                                                 value_type=value_type,
                                                 triu=triu)
                        futures[future] = [cool_type, value_type, chunk_start, output_path]

            temp_zarr_records = []
            for future in as_completed(futures):
                cool_type, value_type, chunk_start, output_path = futures[future]
                try:
                    future.result()
                except Exception as exc:
                    print(f"Cool type {cool_type} value_type {value_type} "
                          f"chunk_start {chunk_start} generated an exception: {exc}")
                    raise exc
                else:
                    temp_zarr_records.append([cool_type, value_type, chunk_start, output_path])

        temp_zarr_records = pd.DataFrame(temp_zarr_records, columns=['cool_type',
                                                                     'value_type',
                                                                     'chunk_start',
                                                                     'output_path'])
        return temp_zarr_records

    @staticmethod
    def _save_single_bin_chunk_worker(ds_paths,
                                      zarr_path,
                                      cool_type,
                                      bin1_slice,
                                      bin2_slice,
                                      sample_id_slice,
                                      value_idx):
        ds = xr.open_mfdataset(ds_paths, concat_dim='sample_id', combine='nested', engine='zarr')
        zarr_da = zarr.open(zarr_path, mode="r+")
        # ds only contains one value type
        _data = ds[cool_type][bin1_slice, bin2_slice, sample_id_slice, 0].load()
        if (_data != 0).sum() == 0:
            # TODO this can be improved, as many bins are empty, but _data.load() still takes significant time
            return

        zarr_da[bin1_slice, bin2_slice, sample_id_slice, value_idx] = _data
        return

    def _save_single_bin_chunk(self, zarr_path, path_array, cool_type):
        with ProcessPoolExecutor(self.n_cpu) as executor:
            futures = {}
            for value_idx, (_, paths) in enumerate(path_array.items()):
                ds = xr.open_mfdataset(paths, concat_dim='sample_id', combine='nested', engine='zarr')
                bin1_idx = ds.get_index('bin1')
                bin2_idx = ds.get_index('bin2')
                sample_id_idx = ds.get_index('sample_id')
                for bin1_start in range(0, bin1_idx.size, self.bin_chunk_size):
                    for bin2_start in range(0, bin2_idx.size, self.bin_chunk_size):
                        for sample_id_start in range(0, sample_id_idx.size, self.sample_chunk_size):
                            bin1_slice = slice(bin1_start, bin1_start + self.bin_chunk_size)
                            bin2_slice = slice(bin2_start, bin2_start + self.bin_chunk_size)
                            sample_id_slice = slice(sample_id_start, sample_id_start + self.sample_chunk_size)
                            future = executor.submit(self._save_single_bin_chunk_worker,
                                                     ds_paths=paths.tolist(),
                                                     zarr_path=zarr_path,
                                                     cool_type=cool_type,
                                                     bin1_slice=bin1_slice,
                                                     bin2_slice=bin2_slice,
                                                     sample_id_slice=sample_id_slice,
                                                     value_idx=value_idx)
                            futures[future] = (bin1_start, bin2_start, sample_id_start)

            for future in as_completed(futures):
                bin1_start, bin2_start, sample_id_start = futures[future]
                try:
                    future.result()
                except Exception as e:
                    print(f'Got error when saving bin1 {bin1_start} to {bin1_start + self.bin_chunk_size} '
                          f'bin2 {bin2_start} to {bin2_start + self.bin_chunk_size} '
                          f'sample_id {sample_id_start} to {sample_id_start + self.sample_chunk_size}')
                    raise e
        return

    def _save_small_sample_chunks_ds_to_final_zarr(self, temp_zarr_records):
        for cool_type, path_df in temp_zarr_records.groupby('cool_type'):
            # make sure the temp zarr paths are in the correct order
            path_array = path_df.pivot(
                index='chunk_start', columns='value_type', values='output_path'
            ).sort_index()[self.value_types[cool_type]]

            zarr_path = f'{self.path}/{cool_type}/'
            self._save_single_bin_chunk(zarr_path=zarr_path,
                                        path_array=path_array,
                                        cool_type=cool_type)
        return

    def _write_data(self):
        # config blosc compressor when run multi-processing
        # see zarr doc here: https://zarr.readthedocs.io/en/stable/tutorial.html#configuring-blosc
        if self.n_cpu > 1:
            from numcodecs import blosc
            blosc.use_threads = False

        final_success_flag = self.log_dir_path / "_WRITE_SUCCESS"
        if final_success_flag.exists():
            print("Data already written. Skipping...")

        # first, convert cool to temp zarr,
        # where the whole matrix is loaded into memory,
        # and sample_id chunk size is small (to limit memory usage)
        temp_zarr_records = self._save_small_sample_chunks()

        # second, write temp zarr to final zarr,
        # where the matrix is chunked into small chunks (self.bin_chunk_size, to limit memory usage),
        # and sample_id chunk size is final (self.sample_chunk_size)
        self._save_small_sample_chunks_ds_to_final_zarr(temp_zarr_records=temp_zarr_records)

        # remove temp zarr files
        for _, row in temp_zarr_records.iterrows():
            shutil.rmtree(row['output_path'])

        # touch success file
        final_success_flag.touch()
        return

    def execute(self):
        """Execute the pipeline."""
        self._init_zarr()
        self._write_data()
        return


def generate_cool_ds(output_dir,
                     cool_table_path,
                     value_types,
                     chrom_sizes_path,
                     trans_matrix=False,
                     mode='w',
                     cooler_bin_size=10000,
                     bin_chunk_size=510,
                     sample_chunk_size=50,
                     data_dtype='float32',
                     cpu=1):
    """
    Generate a CoolDS zarr dataset from cool files.

    Parameters
    ----------
    output_dir :
        The output directory.
    cool_table_path :
        Path to the cool table with four columns: sample, value_type, path, cool_type
    value_types :
        Dict of cool types and their value types.
    chrom_sizes_path :
        Path to the chrom sizes file.
    trans_matrix :
        Whether generate trans-contacts (chrom1 != chrom2) matrix
    mode :
        Mode to open the zarr.
    cooler_bin_size :
        Cooler bin size.
    bin_chunk_size :
        Chunk size of the bin1 and bin2 dimensions.
    sample_chunk_size :
        Chunk size of the sample dimension.
    data_dtype :
        Data type of the matrix.
    cpu :
        Number of CPUs to use.
    """
    # change compressor
    # according to https://zarr.readthedocs.io/en/stable/tutorial.html#compressors
    # and my test, zstd has 1.6x compression ratio compared to default blosc lz4
    compressor = Blosc(cname='zstd', clevel=COMPRESSOR_C_LEVEL, shuffle=Blosc.SHUFFLE)
    zarr.storage.default_compressor = compressor

    output_dir = pathlib.Path(output_dir).absolute().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    chrom_sizes = read_chromsizes(chrom_sizes_path)

    for chrom1 in chrom_sizes.index:
        for chrom2 in chrom_sizes.index:
            if not trans_matrix:
                if chrom1 != chrom2:
                    continue
            if chrom1 == chrom2:
                path = output_dir / chrom1
            else:
                path = output_dir / f'{chrom1}-{chrom2}'

            CoolDSSingleMatrixWriter(path,
                                     cool_table_path=cool_table_path,
                                     value_types=value_types,
                                     chrom_sizes_path=chrom_sizes_path,
                                     chrom1=chrom1,
                                     chrom2=chrom2,
                                     mode=mode,
                                     cooler_bin_size=cooler_bin_size,
                                     bin_chunk_size=bin_chunk_size,
                                     sample_chunk_size=sample_chunk_size,
                                     data_dtype=data_dtype,
                                     cpu=cpu)
    return
