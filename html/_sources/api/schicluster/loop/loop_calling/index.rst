:py:mod:`schicluster.loop.loop_calling`
=======================================

.. py:module:: schicluster.loop.loop_calling


Module Contents
---------------

.. py:function:: fetch_chrom(cool, chrom) -> numpy.array


.. py:function:: select_loop_candidates(cool_e, min_dist, max_dist, resolution, chrom)

   Select loop candidate pixel to perform t test


.. py:function:: paired_t_test(cool_t, cool_t2, chrom, loop, n_cells)

   Paired t test per pixel


.. py:function:: scan_kernel(E, kernel, loop)

   Scan loop surrounding background kernel


.. py:function:: loop_background(E, pad, gap, loop)

   Calculate loop surrounding background level


.. py:function:: call_loop_single_chrom(group_prefix, chrom, resolution=10000, min_dist=50000, max_dist=10000000, pad=5, gap=2)

   calculate t test and loop background for one chromosome


.. py:function:: filter_by_background(data, thres_bl, thres_donut, thres_h, thres_v, resolution)


.. py:function:: find_summit(loop, res, dist_thres)


.. py:function:: call_loops(group_prefix, resolution, output_prefix, thres_bl=1.33, thres_donut=1.33, thres_h=1.2, thres_v=1.2, fdr_thres=0.1, dist_thres=20000, size_thres=1)


.. py:function:: filter_loops(total_loops, output_prefix, fdr_thres, resolution, dist_thres, size_thres)


