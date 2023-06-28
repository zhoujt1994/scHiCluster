:py:mod:`schicluster.loop.merge_cell_to_group`
==============================================

.. py:module:: schicluster.loop.merge_cell_to_group


Module Contents
---------------

.. py:function:: merge_cells_for_single_chromosome(output_dir, output_prefix, merge_type='E')

   Merge cell's E and T matrix to group matrices, sum only, not normalized by n_cells yet:
   E: Matrix normalized by global diagonal backgrounds, calculated from loop_bkg
   E2: E^2 of T, used to calculate global p values
   T: Matrix normalized by global diagonal and local backgrounds, then minus E (T is the delta matrix),
   calculated from loop_bkg
   T2: T^2 of T, used to calculate t test p values


.. py:function:: read_single_cool_chrom(cool_path, chrom, chrom2=None)


.. py:function:: chrom_sum_iterator(chunk_dirs, chrom_sizes, chrom_offset, matrix_type, total_cells)


.. py:function:: save_single_matrix_type(cooler_path, bins_df, chunk_dirs, chrom_sizes, chrom_offset, matrix_type, total_cells)


.. py:function:: merge_group_chunks_to_group_cools(chrom_size_path, resolution, group, output_dir, matrix_types=('E', 'E2', 'T', 'T2', 'Q', 'Q2'))

   Sum all the chunk sum cool files,
   and finally divide the total number of cells to
   get a group cell number normalized cool file in the end.


