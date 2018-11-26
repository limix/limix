*************
API reference
*************

I/O module
==========

.. autosummary::
    :toctree: _generated/

    limix.io.bgen.read
    limix.io.bimbam.read_phenotype
    limix.io.bimbam.see_phenotype
    limix.io.csv.read
    limix.io.csv.see
    limix.io.gen.read
    limix.io.hdf5.fetch
    limix.io.hdf5.fetcher
    limix.io.hdf5.read_limix
    limix.io.hdf5.see
    limix.io.npy.read
    limix.io.npy.save
    limix.io.npy.see
    limix.io.plink.read
    limix.io.plink.see_bed
    limix.io.plink.see_kinship
    limix.io.plink.fetch_dosage

Quality control
===============

.. autosummary::
    :toctree: _generated/

    limix.qc.boxcox
    limix.qc.compute_maf
    limix.qc.count_missingness
    limix.qc.indep_pairwise
    limix.qc.mean_impute
    limix.qc.mean_standardize
    limix.qc.normalise_covariance
    limix.qc.quantile_gaussianize
    limix.qc.regress_out
    limix.qc.remove_dependent_cols
    limix.qc.unique_variants

Statistics
==========

.. autosummary::
    :toctree: _generated/

    limix.stats.allele_expectation
    limix.stats.allele_frequency
    limix.stats.Chi2Mixture
    limix.stats.compute_dosage
    limix.stats.confusion_matrix
    limix.stats.effsizes_se
    limix.stats.empirical_pvalues
    limix.stats.linear_kinship
    limix.stats.lrt_pvalues
    limix.stats.multipletests
    limix.stats.pca

Heritability estimation
=======================

.. autosummary::
    :toctree: _generated/

    limix.her.estimate

Quantitative trait loci
=======================

.. autosummary::
    :toctree: _generated/

    limix.qtl.st_scan
    limix.qtl.QTLModel

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

Generalised Linear Mixed Models
===============================

.. autosummary::
    :toctree: _generated/

    limix.glmm.GLMMComposer.covariance_matrices
    limix.glmm.GLMMComposer.decomp
    limix.glmm.GLMMComposer.fit
    limix.glmm.GLMMComposer.fixed_effects
    limix.glmm.GLMMComposer.likname
    limix.glmm.GLMMComposer.lml
    limix.glmm.GLMMComposer.y

Shell utilities
===============

.. autosummary::
    :toctree: _generated/

    limix.sh.filehash
    limix.sh.download
    limix.sh.extract
    limix.sh.remove
