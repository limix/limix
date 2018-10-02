**********************
Command line interface
**********************

.. plot::
    :context:

    from matplotlib import pyplot as plt
    from limix._cli import cli
    from click.testing import CliRunner
    def call(args):
        CliRunner().invoke(cli, args)
    
    call(["download", "http://rest.s3for.me/limix/small_example.grm.raw.bz2"])
    call(["extract", "small_example.grm.raw.bz2"])
    call(["download", "http://rest.s3for.me/limix/dali.jpg.bz2"])
    call(["extract", "dali.jpg.bz2"])

.. command-output:: limix download http://rest.s3for.me/limix/plink_example.tar.gz && \
      limix extract plink_example.tar.gz &&\
      limix download http://rest.s3for.me/limix/small_example.hdf5 &&\
      limix download http://rest.s3for.me/limix/small_example.csv.bz2 &&\
      limix extract small_example.csv.bz2
    :shell:
    :cwd: _build

Introduction
============

Limix now provides a couple of its functionalities via command line.

.. command-output:: limix --help

You can quickly explore common file types used in genetics, as examples given bellow
will demonstrate.

Kinship
=======

Heatmap representing a plink_ kinship matrix.
Setup::

    limix download http://rest.s3for.me/limix/small_example.grm.raw.bz2
    limix extract small_example.grm.raw.bz2

Command::

    limix see small_example.grm.raw

.. plot::
    :context:
    :include-source: False
    
    call(["see", "small_example.grm.raw"])
    plt.show()


Plink BED format
================

A preview of Plink files in BED format can be done via

.. command-output:: limix -q see plink_example
    :cwd: _build

HDF5
====

The following command shows the hierarchy of a HDF5 file:

.. command-output:: limix see small_example.hdf5
    :cwd: _build

CSV
===

CSV files have their delimiter automatically detected and a preview can be
shown as

.. command-output:: limix see small_example.csv
    :cwd: _build

Image
=====

Finally, an image can be seen via

.. command-output:: limix -q see dali.jpg
    :cwd: _build

.. plot::
    :include-source: False
    :context: close-figs

    >>> call(["see", "dali.jpg"])
    >>> plt.show()

.. _plink: https://www.cog-genomics.org/plink2
.. _hdf5: https://support.hdfgroup.org/HDF5/