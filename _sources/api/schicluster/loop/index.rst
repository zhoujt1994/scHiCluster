:py:mod:`schicluster.loop`
==========================

.. py:module:: schicluster.loop

.. autoapi-nested-parse::

   loop_bkg.py corresponding to draft/loop_bkg_cell, which calculates E and T from cell imputed matrix
   merge_cell_to_group.py corresponding to loop_sumcell_chr,
   which aggregates the cell E, T to cluster level E, E2, T, T2, O,



Submodules
----------
.. toctree::
   :titlesonly:
   :maxdepth: 1

   compare_loop/index.rst
   loop_bkg/index.rst
   loop_calling/index.rst
   merge_cell_to_group/index.rst
   merge_group/index.rst
   merge_raw_matrix/index.rst
   shuffle_fdr/index.rst
   snakemake/index.rst


Package Contents
----------------

.. py:function:: call_loop(cell_table_path, output_dir, chrom_size_path, shuffle=True, chunk_size=200, dist=5050000, cap=5, pad=5, gap=2, resolution=10000, min_cutoff=1e-06, keep_cell_matrix=False, cpu_per_job=10, log_e=True, raw_resolution_str=None, downsample_shuffle=None, black_list_path=None, fdr_pad=7, fdr_min_dist=5, fdr_max_dist=500, fdr_thres=0.1, dist_thres=20000, size_thres=1, cleanup=True)


