:py:mod:`schicluster.loop.loop_bkg`
===================================

.. py:module:: schicluster.loop.loop_bkg


Module Contents
---------------

.. py:function:: calc_diag_stats(E, n_dims)

   Calculate cutoff, average, std, count of non-zero pixels of each diagonals of the E


.. py:function:: calculate_chrom_background_normalization(cell_url, chrom, resolution, output_prefix, dist=5050000, cap=5, pad=5, gap=2, min_cutoff=1e-06, log_e=False, shuffle=False)

   Compute the background for each chromosome in each cell

   :param cell_url:
   :param chrom:
   :param resolution:
   :param output_prefix:
   :param dist:
   :param cap:
   :param pad:
   :param gap:
   :param min_cutoff:
   :param log_e:
   :param shuffle:

   :returns: * *E is the global diagonal normalized matrix*
             * *T is the local background normalized version of E*


