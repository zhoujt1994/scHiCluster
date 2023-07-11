:py:mod:`schicluster.compartment.call_compartment`
==================================================

.. py:module:: schicluster.compartment.call_compartment


Module Contents
---------------

.. py:function:: get_cpg_profile(fasta_path, hdf_output_path, cell_url=None, chrom_size_path=None, resolution=100000)


.. py:function:: compartment_strength(matrix, comp, cpg_ratio)


.. py:function:: single_chrom_compartment(matrix, cpg_ratio, calc_strength=False)


.. py:function:: single_cell_compartment(cell_url, cpg_profile, calc_strength, output_prefix, mode, resolution, chrom_sizes, chrom1, pos1, chrom2, pos2)


.. py:function:: aggregate_compartment(cell_table, temp_dir, bins, output_path, calc_strength)


.. py:function:: multiple_cell_compartment(cell_table_path, output_prefix, cpg_profile_path, cpu=10, calc_strength=False, mode='cool', chrom_size_path=None, resolution=100000, chrom1=1, pos1=2, chrom2=5, pos2=6)


