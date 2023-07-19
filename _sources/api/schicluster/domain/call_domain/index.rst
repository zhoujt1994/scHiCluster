:py:mod:`schicluster.domain.call_domain`
========================================

.. py:module:: schicluster.domain.call_domain


Module Contents
---------------

.. py:data:: PACKAGE_DIR

   

.. py:function:: install_r_package(name)


.. py:function:: domain_df_to_boundary(cool, total_results, resolution)


.. py:function:: single_chrom_calculate_insulation_score(matrix, window_size=10)


.. py:function:: call_domain_and_insulation(cell_url, output_prefix, resolution=25000, window_size=10)


.. py:function:: aggregate_boundary(cell_table, temp_dir, bins, output_path)


.. py:function:: aggregate_insulation(cell_table, temp_dir, bins, output_path)


.. py:function:: multiple_call_domain_and_insulation(cell_table_path, output_prefix, resolution=25000, window_size=10, cpu=10)


