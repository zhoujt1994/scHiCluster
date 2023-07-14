:py:mod:`schicluster.loop.merge_group`
======================================

.. py:module:: schicluster.loop.merge_group


Module Contents
---------------

.. py:function:: chrom_ave_iterator(chunk_dirs, chrom_sizes, chrom_offset, matrix_type, total_cells)


.. py:function:: save_single_matrix_type(cooler_path, bins_df, chunk_dirs, chrom_sizes, chrom_offset, matrix_type, total_cells)


.. py:function:: merge_group_to_bigger_group_cools(chrom_size_path, resolution, group, output_dir, group_list, shuffle, matrix_types=('E', 'E2', 'T', 'T2', 'Q', 'Q2'))

   Sum all the chunk sum cool files,
   and finally divide the total number of cells to
   get a group cell number normalized cool file in the end.


