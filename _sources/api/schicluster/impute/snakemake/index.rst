:py:mod:`schicluster.impute.snakemake`
======================================

.. py:module:: schicluster.impute.snakemake


Module Contents
---------------

.. py:data:: PACKAGE_DIR

   

.. py:function:: prepare_impute(input_scool, output_dir, chrom_size_path, output_dist, window_size, step_size, resolution, batch_size=100, logscale=False, pad=1, std=1, rp=0.5, tol=0.01, min_cutoff=1e-05, cpu_per_job=10)

   prepare snakemake files and directory structure for cell contacts imputation


