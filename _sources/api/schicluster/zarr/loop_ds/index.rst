:py:mod:`schicluster.zarr.loop_ds`
==================================

.. py:module:: schicluster.zarr.loop_ds


Module Contents
---------------

.. py:function:: load_cool_ds_chrom(paths, chrom)


.. py:function:: init_empty_loop_array(cool_ds, loop_mask, value_types, da_name, output_path, chrom, loop_chunk_size=50000)


.. py:function:: save_sample_chunk(output_path, cool_ds_paths, loop_position_ds_path, sample_chunk, sample_start, chrom, da_name, value_types)


.. py:function:: create_loop_ds(cool_ds_paths, loop_position_ds_path, output_path, da_name, chroms, value_types, min_loop_count=1, loop_chunk_size=50000)


