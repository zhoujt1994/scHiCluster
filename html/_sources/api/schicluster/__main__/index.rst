:py:mod:`schicluster.__main__`
==============================

.. py:module:: schicluster.__main__

.. autoapi-nested-parse::

   CLI defined here

   When adding new function:
   1. add a func_register_subparser function to register the subparser
   2. add a condition in main func about this new func name, import the real func as func in main



Module Contents
---------------

.. py:data:: log

   

.. py:data:: DESCRIPTION
   :value: Multiline-String

    .. raw:: html

        <details><summary>Show Value</summary>

    .. code-block:: python

        """
        scHiCluster is a toolkit for single-cell HiC data preprocessing, imputation, and clustering analysis.
        
        Current Tool List in scHiCluster:
        
        """

    .. raw:: html

        </details>

   

.. py:data:: EPILOG
   :value: ''

   

.. py:class:: NiceFormatter(fmt=None, datefmt=None, style='%')


   Bases: :py:obj:`logging.Formatter`

   From Cutadapt https://github.com/marcelm/cutadapt
   Do not prefix "INFO:" to info-level log messages (but do it for all other
   levels).
   Based on http://stackoverflow.com/a/9218261/715090 .

   .. py:method:: format(record)

      Format the specified record as text.

      The record's attribute dictionary is used as the operand to a
      string formatting operation which yields the returned string.
      Before formatting the dictionary, a couple of preparatory steps
      are carried out. The message attribute of the record is computed
      using LogRecord.getMessage(). If the formatting string uses the
      time (as determined by a call to usesTime(), formatTime() is
      called to format the event time. If there is exception information,
      it is formatted using formatException() and appended to the message.



.. py:function:: validate_environment()


.. py:function:: setup_logging(stdout=False, quiet=False, debug=False)

   From Cutadapt https://github.com/marcelm/cutadapt
   Attach handler to the global logger object


.. py:function:: _str_to_bool(v: str) -> bool


.. py:function:: comp_cpg_cell_register_subparser(subparser)


.. py:function:: comp_concatcell_chr_register_subparser(subparser)


.. py:function:: domain_insulation_cell_register_subparser(subparser)


.. py:function:: domain_concatcell_chr_register_subparser(subparser)


.. py:function:: embed_concatcell_chr_register_subparser(subparser)


.. py:function:: embed_mergechr_register_subparser(subparser)


.. py:function:: generatematrix_cell_register_subparser(subparser)


.. py:function:: impute_cell_register_subparser(subparser)


.. py:function:: loop_bkg_cell_register_subparser(subparser)


.. py:function:: loop_sumcell_chr_register_subparser(subparser)


.. py:function:: loop_mergechr_register_subparser(subparser)


.. py:function:: generate_scool_register_subparser(subparser)


.. py:function:: prepare_imputation_register_subparser(subparser)


.. py:function:: call_domain_register_subparser(subparser)


.. py:function:: call_compartment_register_subparser(subparser)


.. py:function:: cpg_ratio_register_subparser(subparser)


.. py:function:: embedding_register_subparser(subparser)


.. py:function:: gene_score_register_subparser(subparser)


.. py:function:: merge_cell_raw_register_subparser(subparser)


.. py:function:: filter_contacts_register_subparser(subparser)


.. py:function:: contact_distance_register_subparser(subparser)


.. py:function:: main()


