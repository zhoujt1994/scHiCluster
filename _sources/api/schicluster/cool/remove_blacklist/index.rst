:py:mod:`schicluster.cool.remove_blacklist`
===========================================

.. py:module:: schicluster.cool.remove_blacklist


Module Contents
---------------

.. py:function:: prepare_2d_blacklist_dict(blacklist_bedpe, resolution=10000)


.. py:function:: _is_2d_blacklist(row, blacklist_2d)


.. py:function:: filter_contacts(contact_path, chrom_size_path=None, blacklist_1d_path=None, blacklist_2d_path=None, output_path=None, remove_duplicates=True, resolution_2d=10000, min_pos_dist=0, chrom1=1, pos1=2, chrom2=5, pos2=6)


.. py:function:: filter_contacts_wrapper(contact_table=None, output_dir=None, chrom_size_path=None, blacklist_1d_path=None, blacklist_2d_path=None, remove_duplicates=True, resolution_2d=10000, chr1=1, pos1=2, chr2=5, pos2=6, min_pos_dist=0, cpu=20)


