:py:mod:`schicluster.cool`
==========================

.. py:module:: schicluster.cool


Submodules
----------
.. toctree::
   :titlesonly:
   :maxdepth: 1

   merge/index.rst
   remove_blacklist/index.rst
   scool/index.rst
   utilities/index.rst


Package Contents
----------------

.. py:function:: generate_scool(contacts_table, output_prefix, chrom_size_path, resolutions, blacklist_1d_path=None, blacklist_2d_path=None, blacklist_resolution=10000, remove_duplicates=True, chr1=1, chr2=5, pos1=2, pos2=6, min_pos_dist=2500, cpu=1, batch_n=50)

   Generate single-resolution cool files from single-cell contact files recorded in contacts_table

   :param contacts_table: tab-separated table containing tow columns, 1) cell id, 2) cell contact file path (juicer-pre format)
                          No header
   :param output_prefix: Output prefix of the cool files. Output path will be {output_prefix}.{resolution_str}.cool
   :param chrom_size_path: Path to the chromosome size file, this file should be in UCSC chrom.sizes format. We will use this file as
                           the reference to create matrix. It is recommended to remove small contigs or chrM from this file.
   :param resolutions: Resolutions to generate the matrix. Each resolution will be stored in a separate file.
   :param blacklist_1d_path: Path to blacklist region BED file, such as ENCODE blacklist.
                             Either side of the contact overlapping with a blacklist region will be removed.
   :param blacklist_2d_path: Path to blacklist region pair BEDPE file.
                             Both side of the contact overlapping with the same blacklist region pair will be removed.
   :param blacklist_resolution: Resolution in bps when consider the 2D blacklist region pairs.
   :param remove_duplicates: If true, will remove duplicated contacts based on [chr1, pos1, chr2, pos2] values
   :param chr1: 0 based index of chr1 column.
   :param chr2: 0 based index of chr2 column.
   :param pos1: 0 based index of pos1 column.
   :param pos2: 0 based index of pos2 column.
   :param min_pos_dist: Minimum distance for a fragment to be considered.
   :param cpu: number of cpus to parallel.
   :param batch_n: number of cells to deal with in each cpu process.


.. py:function:: get_chrom_offsets(bins_df)


.. py:function:: write_coo(path, matrix, chunk_size=5000000)

   Write chromosome contacts as chunked pixel df (cooler input)


.. py:function:: chrom_iterator(input_dir, chrom_order, chrom_offset, chrom_wildcard='{chrom}.hdf', csr=False)


.. py:function:: aggregate_chromosomes(chrom_size_path, resolution, input_dir, output_path, chrom_wildcard='{chrom}.hdf', csr=False)


.. py:function:: cell_chunk(cell_url, chrom_sizes, chunk=50000000)


.. py:function:: aggregate_cells(output_path, cell_dir, chrom_size_path, resolution)


