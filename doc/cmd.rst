**********************
Command line interface
**********************

Limix now provides many of its functionalities via command line.

Quickly explore files
^^^^^^^^^^^^^^^^^^^^^

Show Kinship matrix from a PLINK file:

.. code-block:: bash

    $ limix see some-plink.grm.raw
    # a heatmap should pop-up

The above command will plot a heatmap representing the Kinship matrix stored
in that file.

The following command shows the hierarchy of a HDF5 file:

.. code-block:: bash

    $ limix see 1000G_stage1.hdf5
    /
      +--geauvadis_variants
         +--chr1
         |  +--col_header
         |  |  +--chrom [float64, (951457,)]
         |  |  +--pos [float64, (951457,)]
         |  |  +--rs_ID [|S22, (951457,)]
         |  |  +--snp_type [|S5, (951457,)]
         |  +--matrix [float64, (465, 951457)]
         |  +--row_header
         |     +--sample_ID [|S7, (930,)]

         +--chr10
         |  +--col_header
         |  |  +--chrom [float64, (613471,)]
         |  |  +--pos [float64, (613471,)]
         |  |  +--rs_ID [|S22, (613471,)]
         |  |  +--snp_type [|S5, (613471,)]
         |  +--matrix [float64, (465, 613471)]
         |  +--row_header
         |     +--sample_ID [|S7, (930,)]
         +--chr11
         |  +--col_header
         |  |  +--chrom [float64, (606080,)]
         |  |  +--pos [float64, (606080,)]
         |  |  +--rs_ID [|S22, (606080,)]
         |  |  +--snp_type [|S5, (606080,)]
         |  +--matrix [float64, (465, 606080)]
         |  +--row_header
         |     +--sample_ID [|S7, (930,)]
