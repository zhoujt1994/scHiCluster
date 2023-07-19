:py:mod:`schicluster.loop.merge_cell_to_group`
==============================================

.. py:module:: schicluster.loop.merge_cell_to_group


Module Contents
---------------

.. py:function:: merge_cells_for_single_chromosome(output_dir, output_prefix, merge_type='E')

   Merge cell's E and T matrix to group matrices, sum only, not normalized by n_cells yet
   E: Matrix normalized by global diagonal backgrounds, calculated from loop_bkg
   E2: Sum of square of E, used to calculate global t statistics
   T: Matrix normalized by global diagonal and local backgrounds, then minus E (T is the delta matrix),
   calculated from loop_bkg
   T2: Sum of square of T, used to calculate local t statistics


.. py:function:: read_single_cool_chrom(cool_path, chrom, chrom2=None)


.. py:function:: chrom_sum_iterator(input_cool_list, chrom_sizes, chrom_offset, total_cells)


.. py:function:: save_single_matrix_type(input_cool_list, output_cool, bins_df, chrom_sizes, chrom_offset, total_cells)


.. py:function:: merge_cool(input_cool_tsv_file, output_cool)


.. py:function:: merge_group_chunks_to_group_cools(chrom_size_path, resolution, group, output_dir, matrix_types=('E', 'E2', 'T', 'T2', 'Q', 'Q2'))


.. py:function:: merge_group_to_bigger_group_cools(chrom_size_path, resolution, group, output_dir, group_list, shuffle, matrix_types=('E', 'E2', 'T', 'T2', 'Q', 'Q2'))

   Sum all the group average cool files,
   and finally divide the total number of cells to
   get a group cell number normalized cool file in the end.


