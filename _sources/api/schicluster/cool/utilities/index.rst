:py:mod:`schicluster.cool.utilities`
====================================

.. py:module:: schicluster.cool.utilities


Module Contents
---------------

.. py:function:: get_chrom_offsets(bins_df)


.. py:function:: write_coo(path, matrix, chunk_size=5000000)

   Write chromosome contacts as chunked pixel df (cooler input)


.. py:function:: chrom_iterator(input_dir, chrom_order, chrom_offset, chrom_wildcard='{chrom}.hdf', csr=False)


.. py:function:: aggregate_chromosomes(chrom_size_path, resolution, input_dir, output_path, chrom_wildcard='{chrom}.hdf', csr=False)


.. py:function:: cell_chunk(cell_url, chrom_sizes, chunk=50000000)


.. py:function:: aggregate_cells(output_path, cell_dir, chrom_size_path, resolution)


