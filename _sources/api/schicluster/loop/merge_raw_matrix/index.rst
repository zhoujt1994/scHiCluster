:py:mod:`schicluster.loop.merge_raw_matrix`
===========================================

.. py:module:: schicluster.loop.merge_raw_matrix


Module Contents
---------------

.. py:function:: _chrom_sum_iterator(cell_urls, chrom_sizes, chrom_offset, add_trans=False)

   Iterate through the raw matrices and chromosomes of cells.

   :param cell_urls: List of cell urls.
   :param chrom_sizes: Dictionary of chromosome sizes.
   :param chrom_offset: Dictionary of chromosome offsets.
   :param add_trans: If true, will also iterate all the trans combinations (different chromosomes).

   :Yields: *pixel_df* -- Dataframe of pixels. Used by cooler.create_cooler to save to h5 file.


.. py:function:: _save_single_matrix_type(cooler_path, bins_df, cell_urls, chrom_sizes, chrom_offset, add_trans=False)

   Save a single matrix type Cool file from merging multiple cell urls.

   :param cooler_path: Path to the output cool file.
   :param bins_df: Dataframe of bins. Created from chromosome sizes and resolution.
   :param cell_urls: List of cell urls to merge.
   :param chrom_sizes: Dictionary of chromosome sizes.
   :param chrom_offset: Dictionary of chromosome offsets.
   :param add_trans: Whether add trans matrix also.


.. py:function:: make_raw_matrix_cell_table(cell_table, resolution_str='10K')


.. py:function:: merge_raw_scool_by_cluster(chrom_size_path, resolution, cell_table_path, output_dir, add_trans=False, cpu=1)

   Sum the raw matrix of cells, no normalization.

   :param chrom_size_path: Path to the chrom size file. This file is used to determine chromosome names and bins.
   :param resolution: Resolution of the raw matrix.
   :param cell_table_path: Path to the cell table.
                           This table should contain three columns: cell_id, cell_url, cell_group; no Header.
                           The cell_id is the id of the cell in the raw matrix.
                           The cell_url is the path to the raw matrix.
                           The cell_group is the group of the cells.
   :param output_dir: Path to the output directory. Group cool files will be named as "<output_dir>/<cell_group>.cool".
   :param add_trans: Whether add trans matrix also.
   :param cpu: Number of CPUs to use.


