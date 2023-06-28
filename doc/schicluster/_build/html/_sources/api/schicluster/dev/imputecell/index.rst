:py:mod:`schicluster.dev.imputecell`
====================================

.. py:module:: schicluster.dev.imputecell


Module Contents
---------------

.. py:data:: hg38dim
   :value: [248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636,...

   

.. py:data:: hg19dim
   :value: [249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022,...

   

.. py:data:: mm10dim
   :value: [195471971, 182113224, 160039680, 156508116, 151834684, 149736546, 145441459, 129401213,...

   

.. py:function:: random_walk_cpu(P, rp, tol, dist, spr)


.. py:function:: impute_cell(indir, outdir, cell, chrom, res, genome, mode, logscale=False, pad=1, std=1, rp=0.5, tol=0.01, rwr_dist=500000000, rwr_sparsity=1, output_dist=500000000, output_format='hdf')


