:py:mod:`schicluster.schicluster`
=================================

.. py:module:: schicluster.schicluster


Module Contents
---------------

.. py:function:: neighbor_ave_gpu(A, pad)


.. py:function:: random_walk_gpu(A, rp)


.. py:function:: impute_gpu(args)


.. py:function:: hicluster_gpu(network, chromsize, nc, res=1000000, pad=1, rp=0.5, prct=20, ndim=20)


.. py:function:: neighbor_ave_cpu(A, pad)


.. py:function:: random_walk_cpu(A, rp)


.. py:function:: impute_cpu(args)


.. py:function:: hicluster_cpu(network, chromsize, nc, res=1000000, pad=1, rp=0.5, prct=20, ndim=20, ncpus=10)


.. py:function:: raw_pca(network, chromsize, nc, res=1000000, ndim=20)


.. py:function:: ds_pca(network, chromsize, nc, res=1000000, ndim=20)


.. py:function:: comp_comp(O, cg)


.. py:function:: compartment(network, chromsize, nc, res=1000000, ndim=20)


.. py:function:: decay(network, chromsize, nc, res=1000000, ndim=20)


.. py:function:: merge_gpu(network, c, res, pad=1, rp=0.5, prct=-1)


.. py:function:: merge_cpu(network, c, res, pad=1, rp=0.5, prct=-1)


.. py:function:: output_topdom(cell, c, Q, res)


.. py:function:: output_sparse(cell, c, Q, res)


.. py:function:: diff_dom(args)


.. py:function:: chrflankpv(args)


.. py:function:: filter_bins(prob, fdr, fdr_cutoff=0.05, diff_cutoff=0, max_cutoff=0, min_cutoff=1)


