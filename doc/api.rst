.. currentmodule:: limix

*************
API reference
*************

Example files
=============

.. autosummary::
    :toctree: _generated/

    example.file_example

Generalised Linear Mixed Models
===============================

.. autosummary::
    :toctree: _generated/

    glmm
    glmm.GLMMComposer

Heritability estimation
=======================

.. autosummary::
    :toctree: _generated/

    limix.her.estimate

I/O module
==========

limix.io.bgen
-------------

.. autosummary::
    :toctree: _generated/

    limix.io.bgen.read

limix.io.bimbam
---------------

.. autosummary::
    :toctree: _generated/

    limix.io.bimbam.read_phenotype
    limix.io.bimbam.see_phenotype

limix.io.csv
------------

.. autosummary::
    :toctree: _generated/

    limix.io.csv.read
    limix.io.csv.see

limix.io.gen
------------

.. autosummary::
    :toctree: _generated/

    limix.io.gen.read

limix.io.hdf5
-------------

.. autosummary::
    :toctree: _generated/

    limix.io.hdf5.fetch
    limix.io.hdf5.fetcher
    limix.io.hdf5.read_limix
    limix.io.hdf5.see

limix.io.npy
------------

.. autosummary::
    :toctree: _generated/

    limix.io.npy.read
    limix.io.npy.save
    limix.io.npy.see

limix.io.plink
--------------

.. autosummary::
    :toctree: _generated/

    limix.io.plink.read
    limix.io.plink.see_bed
    limix.io.plink.see_kinship
    limix.io.plink.fetch_dosage

Plotting & Graphics
===================

.. autosummary::
    :toctree: _generated/

    limix.plot.box_aspect
    limix.plot.ConsensusCurve
    limix.plot.image
    limix.plot.kinship
    limix.plot.load_dataset
    limix.plot.manhattan
    limix.plot.normal
    limix.plot.pca
    limix.plot.power
    limix.plot.qqplot
    limix.plot.image
    limix.plot.get_pyplot
    limix.plot.show


Quanlity control
================

.. autosummary::
    :toctree: _generated/

    limix.qc.compute_maf
    limix.qc.boxcox
    limix.qc.mean_standardize
    limix.qc.quantile_gaussianize
    limix.qc.regress_out
    limix.qc.remove_dependent_cols
    limix.qc.mean_impute
    limix.qc.indep_pairwise
    limix.qc.count_missingness
    limix.qc.compute_maf
    limix.qc.normalise_covariance
    limix.qc.unique_variants

Quantitative trait loci
=======================

.. autosummary::
    :toctree: _generated/

    limix.qtl.scan
    limix.qtl.QTLModel

Shell utilities
===============

.. autosummary::
    :toctree: _generated/

    limix.sh.filehash
    limix.sh.download
    limix.sh.extract
    limix.sh.remove


Statistics
==========

.. autosummary::
    :toctree: _generated/

    limix.stats.pca
    limix.stats.multipletests
    limix.stats.empirical_pvalues
    limix.stats.Chi2Mixture
    limix.stats.linear_kinship
    limix.stats.lrt_pvalues
    limix.stats.effsizes_se
    limix.stats.confusion_matrix
    limix.stats.convert_to_dosage
