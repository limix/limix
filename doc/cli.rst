**********************
Command line interface
**********************

.. plot::
   :context:

   >>> from limix.plot import show
   >>> from limix._cli import cli
   >>> from limix.sh import remove as rm
   >>> from click.testing import CliRunner
   >>> def call(args):
   ...     CliRunner().invoke(cli, args, catch_exceptions=False)

Introduction
============

Limix now provides a couple of its functionalities via command line.

.. command-output:: limix --help

You can quickly explore common file types used in genetics, as examples given bellow
will demonstrate.

Kinship
=======

Heatmap representing a plink_ kinship matrix.

.. command-output:: limix download http://rest.s3for.me/limix/small_example.grm.raw.bz2
   :cwd: _build

.. command-output:: limix extract small_example.grm.raw.bz2
   :cwd: _build

.. command-output:: limix see small_example.grm.raw
   :cwd: _build

.. plot::
   :include-source: False
   :context: close-figs

   >>> call(["download", "http://rest.s3for.me/limix/small_example.grm.raw.bz2"])
   >>> call(["extract", "small_example.grm.raw.bz2"])
   >>> call(["see", "small_example.grm.raw"])
   >>> rm("small_example.grm.raw")
   >>> rm("small_example.grm.raw.bz2")

Plink BED format
================

A preview of Plink files in BED format can be done via

.. command-output:: limix download http://rest.s3for.me/limix/plink_example.tar.gz
   :cwd: _build

.. command-output:: limix extract plink_example.tar.gz
   :cwd: _build

.. command-output:: limix see plink_example
   :cwd: _build

.. plot::
   :include-source: False
   :context: close-figs

   >>> call(["download", "http://rest.s3for.me/limix/plink_example.tar.gz"])
   >>> call(["extract", "plink_example.tar.gz"])
   >>> call(["see", "plink_example"])
   >>> rm("plink_example.bed")
   >>> rm("plink_example.bim")
   >>> rm("plink_example.fam")
   >>> rm("plink_example.tar.gz")

BIMBAM file formats
===================

Phenotype:

.. command-output:: limix download http://rest.s3for.me/limix/ex0/phenotype.gemma
   :cwd: _build

.. command-output:: limix see phenotype.gemma:bimbam-pheno
   :cwd: _build

.. plot::
   :include-source: False
   :context: close-figs

   >>> call(["download", "http://rest.s3for.me/limix/ex0/phenotype.gemma"])
   >>> call(["see", "phenotype.gemma:bimbam-pheno"])
   >>> rm("phenotype.gemma")

HDF5
====

The following command shows the hierarchy of a HDF5 file:

.. command-output:: limix download http://rest.s3for.me/limix/small_example.hdf5
   :cwd: _build

.. command-output:: limix see small_example.hdf5
   :cwd: _build

.. plot::
   :include-source: False
   :context: close-figs

   >>> call(["download", "http://rest.s3for.me/limix/small_example.hdf5"])
   >>> call(["see", "small_example.hdf5"])
   >>> rm("small_example.hdf5")

CSV
===

CSV files have their delimiter automatically detected and a preview can be
shown as


.. command-output:: limix download http://rest.s3for.me/limix/small_example.csv.bz2
   :cwd: _build

.. command-output:: limix extract small_example.csv.bz2
   :cwd: _build

.. command-output:: limix see small_example.csv
   :cwd: _build

.. plot::
   :include-source: False
   :context: close-figs

   >>> call(["download", "http://rest.s3for.me/limix/small_example.csv.bz2"])
   >>> call(["extract", "small_example.csv.bz2"])
   >>> call(["see", "small_example.csv"])
   >>> rm("small_example.csv")
   >>> rm("small_example.csv.bz2")

Image
=====

An image can be seen via

.. command-output:: limix download http://rest.s3for.me/limix/dali.jpg.bz2
   :cwd: _build

.. command-output:: limix extract dali.jpg.bz2
   :cwd: _build

.. command-output:: limix see dali.jpg
   :cwd: _build

.. plot::
   :include-source: False
   :context: close-figs

   >>> call(["download", "http://rest.s3for.me/limix/dali.jpg.bz2"])
   >>> call(["extract", "dali.jpg.bz2"])
   >>> call(["see", "dali.jpg"])
   >>> rm("dali.jpg")
   >>> rm("dali.jpg.bz2")

GWAS
====

.. command-output:: limix qtl --help
   :cwd: _build

.. plot::
   :include-source: False
   :context: close-figs

   >>> call(["qtl", "--help"])

.. _plink: https://www.cog-genomics.org/plink2

