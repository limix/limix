*************
API reference
*************

I/O module
==========

.. autosummary::
    :toctree: api/

    limix.io.bgen.read
    limix.io.bimbam.read_phenotype
    limix.io.csv.read
    limix.io.gen.read
    limix.io.hdf5.read_limix
    limix.io.npy.read
    limix.io.plink.read

Quality control
===============

.. autosummary::
    :toctree: api/

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
    :toctree: api/

    limix.stats.Chi2Mixture
    limix.stats.allele_expectation
    limix.stats.allele_frequency
    limix.stats.compute_dosage
    limix.stats.confusion_matrix
    limix.stats.empirical_pvalues
    limix.stats.linear_kinship
    limix.stats.lrt_pvalues
    limix.stats.multipletests
    limix.stats.pca

Heritability estimation
=======================

.. autosummary::
    :toctree: api/

    limix.her.estimate

Variance decomposition
======================

.. autosummary::
    :toctree: api/

    limix.vardec.VarDec

Quantitative trait loci
=======================

.. autosummary::
    :toctree: api/

    limix.qtl.scan
    limix.qtl.iscan

Plotting & Graphics
===================

.. autosummary::
    :toctree: api/

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

Shell utilities
===============

.. autosummary::
    :toctree: api/

    limix.sh.filehash
    limix.sh.download
    limix.sh.extract
    limix.sh.remove
