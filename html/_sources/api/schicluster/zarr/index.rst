:py:mod:`schicluster.zarr`
==========================

.. py:module:: schicluster.zarr


Submodules
----------
.. toctree::
   :titlesonly:
   :maxdepth: 1

   cool_ds/index.rst
   loop_ds/index.rst


Package Contents
----------------

.. py:class:: CoolDSSingleMatrixWriter(path, cool_table_path, value_types, chrom_sizes_path, chrom1, chrom2=None, mode='w', cooler_bin_size=10000, bin_chunk_size=510, sample_chunk_size=50, data_dtype='float32', cpu=1)


   .. py:method:: _read_cool_table(cool_table_path)


   .. py:method:: _read_chrom_info(chrom_sizes_path)


   .. py:method:: _add_root_attrs()


   .. py:method:: _init_zarr()


   .. py:method:: _save_cool_to_temp_zarr(cool_paths, output_path, cool_type, value_type, triu)


   .. py:method:: _save_small_sample_chunks()


   .. py:method:: _save_single_bin_chunk_worker(ds_paths, zarr_path, cool_type, bin1_slice, bin2_slice, sample_id_slice, value_idx)
      :staticmethod:


   .. py:method:: _save_single_bin_chunk(zarr_path, path_array, cool_type)


   .. py:method:: _save_small_sample_chunks_ds_to_final_zarr(temp_zarr_records)


   .. py:method:: _write_data()


   .. py:method:: execute()

      Execute the pipeline.



.. py:function:: generate_cool_ds(output_dir, cool_table_path, value_types, chrom_sizes_path, trans_matrix=False, mode='w', cooler_bin_size=10000, bin_chunk_size=510, sample_chunk_size=50, data_dtype='float32', cpu=1)

   Generate a CoolDS zarr dataset from cool files.

   :param output_dir: The output directory.
   :param cool_table_path: Path to the cool table with four columns: sample, value_type, path, cool_type
   :param value_types: Dict of cool types and their value types.
   :param chrom_sizes_path: Path to the chrom sizes file.
   :param trans_matrix: Whether generate trans-contacts (chrom1 != chrom2) matrix
   :param mode: Mode to open the zarr.
   :param cooler_bin_size: Cooler bin size.
   :param bin_chunk_size: Chunk size of the bin1 and bin2 dimensions.
   :param sample_chunk_size: Chunk size of the sample dimension.
   :param data_dtype: Data type of the matrix.
   :param cpu: Number of CPUs to use.


