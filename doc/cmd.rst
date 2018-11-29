**********************
Command line interface
**********************

.. plot::
    :context:

    >>> from limix.plot import show
    >>> from limix._cli import cli
    >>> from click.testing import CliRunner
    >>> def call(args):
    ...     CliRunner().invoke(cli, args)
    >>>
    >>> call(["download", "http://rest.s3for.me/limix/small_example.grm.raw.bz2"])
    >>> call(["extract", "small_example.grm.raw.bz2"])
    >>> call(["download", "http://rest.s3for.me/limix/dali.jpg.bz2"])
    >>> call(["extract", "dali.jpg.bz2"])

.. command-output:: limix download -q http://rest.s3for.me/limix/plink_example.tar.gz &\
      limix download -q http://rest.s3for.me/limix/small_example.hdf5 &\
      limix download -q http://rest.s3for.me/limix/small_example.csv.bz2 &\
      limix download -q http://rest.s3for.me/limix/ex0/phenotype.gemma
    :shell:
    :cwd: _build

.. command-output:: limix extract -q plink_example.tar.gz &\
      limix extract -q small_example.csv.bz2
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

    >>> call(["see", "small_example.grm.raw"])
    >>> show()


Plink BED format
================

A preview of Plink files in BED format can be done via

.. command-output:: limix see plink_example
    :cwd: _build

BIMBAM file formats
===================

Phenotype:

.. command-output:: limix see phenotype.gemma:bimbam-pheno
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

An image can be seen via

.. command-output:: limix see -q dali.jpg
    :cwd: _build

.. plot::
    :include-source: False
    :context: close-figs

    >>> call(["see", "dali.jpg"])
    >>> show()

GWAS
====

.. command-output:: limix scan --help
    :cwd: _build

.. _plink: https://www.cog-genomics.org/plink2

