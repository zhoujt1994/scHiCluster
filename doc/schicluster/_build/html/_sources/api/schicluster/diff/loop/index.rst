:py:mod:`schicluster.diff.loop`
===============================

.. py:module:: schicluster.diff.loop


Module Contents
---------------

.. py:function:: one_way_anova(chrom_loop_ds, da_name, value_type, group_n_dim='group_n', group_dim='sample_id')

   Perform one-way ANOVA on a single-chrom loop dataset.

   :param chrom_loop_ds: A single-chrom loop dataset.
   :param da_name: The name of the data array to perform ANOVA on.
   :param value_type: The value type of the data array to perform ANOVA on.
                      both "{value_type}" and "{value_type}2" should be present in the "{da_name}_value_type" dimension.
   :param group_n_dim: The name of the group number variable.
   :param group_dim: The name of the group dimension.

   :rtype: F statistics and P-values of the one-way ANOVA.


.. py:function:: merge_groups(loop_ds, group_map, da_name, group_dim='sample_id', group_n_dim='group_n')

   Merge groups into larger groups in a loop dataset.

   :param loop_ds: A loop dataset.
   :param group_map: A pd.Series mapping from old group names to new group names.
   :param da_name: The name of the data array to merge groups for.
   :param group_dim: The name of the group dimension.
   :param group_n_dim: The name of the group number variable.

   :rtype: A loop dataset with merged groups.


