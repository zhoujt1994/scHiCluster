:py:mod:`schicluster.loop.snakemake`
====================================

.. py:module:: schicluster.loop.snakemake


Module Contents
---------------

.. py:data:: PACKAGE_DIR

   

.. py:data:: GENERATE_MATRIX_CHUNK_TEMPLATE

   

.. py:data:: GENERATE_MATRIX_SCOOL_TEMPLATE

   

.. py:function:: prepare_dir(output_dir, chunk_df, dist, cap, pad, gap, resolution, min_cutoff, chrom_size_path, keep_cell_matrix, log_e_str, shuffle)


.. py:function:: prepare_loop_snakemake(cell_table_path, output_dir, chrom_size_path, chunk_size=100, dist=5050000, cap=5, pad=5, gap=2, resolution=10000, min_cutoff=1e-06, keep_cell_matrix=False, cpu_per_job=10, log_e=True, shuffle=False, raw_resolution_str=None, downsample_shuffle=None)


.. py:function:: check_chunk_dir_finish(output_dir)


.. py:function:: _run_snakemake(output_dir)


.. py:function:: call_loop(cell_table_path, output_dir, chrom_size_path, shuffle=True, chunk_size=200, dist=5050000, cap=5, pad=5, gap=2, resolution=10000, min_cutoff=1e-06, keep_cell_matrix=False, cpu_per_job=10, log_e=True, raw_resolution_str=None, downsample_shuffle=None, black_list_path=None, fdr_pad=7, fdr_min_dist=5, fdr_max_dist=500, fdr_thres=0.1, dist_thres=20000, size_thres=1, cleanup=True)


.. py:function:: merge_loop(group, output_dir, chrom_size_path, shuffle=True, chunk_size=200, dist=5050000, cap=5, pad=5, gap=2, resolution=10000, min_cutoff=1e-06, keep_cell_matrix=False, cpu_per_job=10, log_e=True, raw_resolution_str=None, downsample_shuffle=None, black_list_path=None, fdr_pad=7, fdr_min_dist=5, fdr_max_dist=500, fdr_thres=0.1, dist_thres=20000, size_thres=1, cleanup=True)


